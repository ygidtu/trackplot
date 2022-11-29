#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Generate object for plotting reads like IGV track
"""
import os.path
import sys
from typing import Optional, List, Dict, Union

import numpy as np
import pandas as pd
import pysam
from loguru import logger

from sashimi.base.CoordinateMap import Coordinate
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.Readder import Reader
from sashimi.file.File import File

np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


class Reads(GenomicLoci):

    def __init__(self,
                 chromosome: str,
                 start: str,
                 end: str,
                 strand: str,
                 id: str,
                 exons: list,
                 introns: list,
                 polya_length: float = -1.0,
                 m6a: float = -1.0,
                 features: Optional[list] = None):
        u"""
        Fetch information of each read from bam file, for stacked plot like IGV
        :param chromosome: the chromosome id of the given read
        :param start: the start site of the given read
        :param end: the end site of the given read
        :param strand: the strand of the given read
        :param id: the name of the given read
        :param exons: a list of GenomicLoci obj which contains region of exon
        :param introns: a list of GenomicLoci obj which contains region of intron
        :param polya_length: length of polya
        :param features: current support m6a site and polya length, like below
        """

        super().__init__(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand
        )

        self.exons = self.__collapse_read__(sorted(exons))
        self.introns = self.__collapse_read__(sorted(introns))
        self.polya_length = polya_length
        self.m6a = m6a
        self.features = features
        self.id = id

    @property
    def exon_list(self):
        u"""
        return a nested list which contains exon regions
        1-based coordinate
        :return: a nested list
        """

        exon_nested_lst = []
        for i in self.exons:
            exon_nested_lst.append(([i.start, i.end]))
        return exon_nested_lst

    @property
    def intron_list(self):
        u"""
        return a nested list which contains intronic regions
        1-based coordinate
        :return: a nested list
        """

        intron_nested_lst = []
        for i in self.introns:
            intron_nested_lst.append(
                ([i.start, i.end])
            )
        return intron_nested_lst

    def __len__(self):
        return sum(map(lambda x: x[1] - x[0] + 1, self.exon_list))

    def __str__(self):
        exons_str = []
        for i in self.exons:
            exons_str.append("{}-{}".format(i.start, i.end))

        return "{}:{}-{}:{} {} {}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand,
            self.id,
            "|".join(exons_str)
        )

    def to_dict(self) -> Dict:
        u"""
        return a dict for generate DataFrame object
        :return:
        """
        return {
            "chromosome": self.chromosome,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "polya_length": self.polya_length,
            "m6a": self.m6a
        }

    @staticmethod
    def __collapse_read__(genomic_list: List[GenomicLoci]) -> List[GenomicLoci]:
        u"""
        Collapse the consecutive region, because of the coordinate shift of matplotlib
        :param genomic_list: a list of GenomicLoci
        :return:
        """

        if len(genomic_list) <= 1:
            return genomic_list

        list_return = []
        current_gl = None

        for gl in genomic_list:
            if not current_gl:
                current_gl = gl
            else:
                if current_gl.end + 1 == gl.start:
                    current_gl = GenomicLoci(
                        chromosome=current_gl.chromosome,
                        start=current_gl.start,
                        end=gl.end,
                        strand=current_gl.strand,
                        name=current_gl.name
                    )
                else:
                    list_return.append(current_gl)
                    current_gl = gl

        list_return.append(current_gl)
        return list_return


class ReadSegment(File):

    def __init__(
            self,
            path: str,
            label: str = "",
            meta: Optional[pd.DataFrame] = None,
            region: Optional[GenomicLoci] = None,
            library: str = "fru",
            deletion_ignore: Optional[int] = True,
            del_ratio_ignore: float = .5,
            features: Optional[dict] = None,
            exon_focus: Optional[str] = None,
            is_bed: bool = False
    ):
        u"""
        init a class for store the information for IGV-like plot
        :param path: the path of the input file
        :param label: the label of current track
        :param meta: the meta info of all reads which was used to sort the reads
        :param region: the region for plotting
        :param library: the library category of the input file
        :param deletion_ignore: whether to ignore the deletion region, default: Ture
        :param del_ratio_ignore: the length of deletion must below to del_ratio_ignore * length of alignment
        :param features: default: features={"m6a": "ma","real_strand": "rs","polya": "pa"}
        :param exon_focus: exon to focus, like start1-end1,start2-end2
        :param is_bed: bed file for generating igv-like reads.
        """
        super().__init__(path)

        self.features = None
        assert library in ["frf", "frs", "fru"], \
            "Illegal library name."

        self.library = library
        self.data = []
        self.meta = meta
        self.region = region
        self.label = label
        self.deletion_ignore = deletion_ignore
        self.del_ratio_ignore = del_ratio_ignore
        self.features = features
        self.is_bed = is_bed
        self.exon_focus = set(map(lambda x: x.strip(), exon_focus.split(','))) if exon_focus else exon_focus

    @classmethod
    def create(
            cls,
            path: str,
            label: str = "",
            library: str = "fru",
            deletion_ignore: Optional[int] = True,
            del_ratio_ignore: float = .5,
            features: Optional[dict] = None,
            exon_focus: Optional[str] = None

    ):
        u"""
        Load each reads and its strand, features.
        m6a tag support the genomic site
        polya tag support length of polya. if polya tag is available, then read_strand (rl) is also essential

        :param path:
        :param label:
        :param library:
        :param deletion_ignore: ignore the deletion length
        :param del_ratio_ignore: ignore the deletion length which calculated by mapped length * ratio
        :param features: support m6a and polyA length from bam tag.
        like {"m6a": "ma", "polya": "pa", "real_strand": "rs"}
        :param exon_focus: exon to focus, like start1-end1,start2-end2
        :return:
        """
        # if bam file then check the index file of the given file
        if path.endswith("bam") and not os.path.exists(path + ".bai"):
            pysam.index(path)

        is_bed = False

        # if bed file
        if "bed" in path:
            is_bed = True
            if not os.path.exists(path + ".tbi"):
                try:
                    path = pysam.tabix_index(path, preset="bed", force=True)
                except OSError as e:
                    logger.error(f"Failed to build index for {path}. \n {e} \n"
                                 f"Please sort and index your bed file by "
                                 f"`bedtools sort -i {path} | bgzip > {path}.gz`")
                    sys.exit(0)
                    # path_new = re.sub(".bed.gz$", "", path) + 'sorted.bed.gz'
                    #
                    # Reference.sort_gtf(input_gtf=path,
                    #                    output_gtf=path_new)
                    # path = path_new

        return cls(
            path=path,
            label=label,
            library=library,
            deletion_ignore=deletion_ignore,
            del_ratio_ignore=del_ratio_ignore,
            features=features,
            exon_focus=exon_focus,
            is_bed=is_bed
        )

    def get_index(self):
        u"""
        get a nested list which presents the order of plot
        :return:
        """
        assert self.meta is not None, f"Not found meta information, please `load` first."
        for ind in self.meta.groupby(['y_loci'])['list_index'].apply(list).tolist()[::-1]:
            yield ind

    def set_region(self,
                   chromosome,
                   start,
                   end,
                   strand):
        u"""
        set the plotting region
        :param chromosome: the name of chromosome
        :param start: the start of the plotting region
        :param end: the end of the plotting region
        :param strand: stand of the plotting region
        :return:
        """
        self.region = GenomicLoci(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand
        )

    @staticmethod
    def df_sort(dfs: pd.DataFrame) -> Optional[pd.DataFrame]:
        u"""
        sorting the dataframe to generate plot index,
        copy from jinbu jia
        :param dfs: a pd.DataFrame object
        :return: pd.DataFrame
        """

        y_loci = []
        have_overlap_regions = []
        now_max_x = 0
        height = 1
        y_min = 1
        dfs_lst = []
        for _, df in dfs.groupby('exon_group'):
            for index, row in df.iterrows():
                start = row["start"]
                end = row["end"]

                if start > now_max_x:
                    if y_min != 0:
                        y_min += 1
                    y_max = height
                    now_max_x = end
                    have_overlap_regions = [[1, y_max, now_max_x]]
                else:
                    for d in have_overlap_regions:
                        if d[2] < start:
                            d[2] = start - 1
                    if len(have_overlap_regions) > 1:
                        new_have_overlap_regions = [have_overlap_regions[0]]
                        for d in have_overlap_regions[1:]:
                            if d[2] == new_have_overlap_regions[-1][2]:
                                new_have_overlap_regions[-1][1] = d[1]
                            else:
                                new_have_overlap_regions.append(d)
                        have_overlap_regions = new_have_overlap_regions
                    have_insert = False
                    for (i, d) in enumerate(have_overlap_regions):
                        x1, x2, x3 = d
                        if x3 < start and (x2 - x1 + 1) >= height:
                            have_insert = True
                            y_min = x1
                            y_max = y_min + height - 1
                            if y_max != x3:
                                d[0] = y_max + 1
                                have_overlap_regions.insert(i, [x1, y_max, end])
                            else:
                                d[2] = end
                            break
                    if not have_insert:
                        y_min = have_overlap_regions[-1][1] + 1 if have_overlap_regions else 1
                        y_max = y_min + height - 1
                        have_overlap_regions.append([y_min, y_max, end])
                    if end > now_max_x:
                        now_max_x = end
                y_loci.append(y_min)
            df["y_loci"] = y_loci
            dfs_lst.append(df)
            y_loci = []

        if len(dfs_lst) <= 0:
            logger.error(f"There is no any read segments")
            return None

        return pd.concat(dfs_lst)

    def load_bed(self):
        try:
            for rec in Reader.read_gtf(self.path, region=self.region, bed=True):
                exon_bound = []
                intron_bound = []
                current_start = int(rec[1])
                current_end = int(rec[2])
                if len(rec) > 3:
                    current_id = rec[3]
                else:
                    current_id = "NoID"

                if len(rec) != 12:
                    exon_bound.append(
                        GenomicLoci(
                            chromosome=self.region.chromosome,
                            start=current_start + 1,
                            end=current_end,
                            strand=self.region.strand,
                            name="exon"
                        )
                    )
                else:

                    block_sizes = [int(x) for x in rec[10].split(",") if x]
                    block_starts = [int(x) for x in rec[11].split(",") if x]

                    read = []
                    for i in range(len(block_starts)):
                        exon_bound.append(
                            GenomicLoci(
                                chromosome=self.region.chromosome,
                                start=current_start + 1 + block_starts[i],
                                end=current_start + 1 + block_starts[i] + block_sizes[i] - 1,
                                strand=self.region.strand,
                                name="exon"
                            )
                        )
                    for pre_exon, next_exon in Coordinate.__slide_window__(exon_bound, 2):
                        intron_bound.append(
                            GenomicLoci(
                                chromosome=self.region.chromosome,
                                start=pre_exon.end + 1,
                                end=next_exon.start - 1,
                                strand=self.region.strand,
                                name="intron"
                            )
                        )

                read = Reads(
                    chromosome=self.region.chromosome,
                    start=min(map(lambda x: x.start, exon_bound)),
                    end=max(map(lambda x: x.end, exon_bound)),
                    strand=self.region.strand,
                    id=current_id,
                    exons=exon_bound,
                    introns=intron_bound,
                    polya_length=-1,
                    m6a=-1,
                    features=[]
                )
                if read.start < self.region.start or read.end > self.region.end:
                    continue

                self.data.append(read)

        except IOError as err:
            logger.error('There is no .bed file at {0}'.format(self.path))
            logger.error(err)
        except ValueError as err:
            logger.error(self.path)
            logger.error(err)

    def load_bam(self):
        try:
            for read, _ in Reader.read_bam(self.path, self.region):
                if read.reference_start < self.region.start or read.reference_end > self.region.end:
                    continue

                if not self.deletion_ignore:
                    current_ignore_num = min(
                        [self.deletion_ignore, read.query_alignment_length * self.del_ratio_ignore])
                else:
                    current_ignore_num = np.Inf

                exon_bound = []
                intron_bound = []

                start = read.reference_start + 1
                for c, l in read.cigar:
                    if c == 0:  # for match
                        exon_bound.append(
                            GenomicLoci(
                                chromosome=self.region.chromosome,
                                start=start,
                                end=start + l - 1,
                                strand=self.region.strand,
                                name="exon"
                            )
                        )
                        start += l

                    elif c == 1:  # for insert
                        continue

                    elif c == 2:  # for del
                        if l <= current_ignore_num:
                            exon_bound.append(
                                GenomicLoci(
                                    chromosome=self.region.chromosome,
                                    start=start,
                                    end=start + l - 1,
                                    strand=self.region.strand,
                                    name="exon"
                                )
                            )
                        start += l

                    elif c == 3:  # for intron
                        intron_bound.append(
                            GenomicLoci(
                                chromosome=self.region.chromosome,
                                start=start,
                                end=start + l - 1,
                                strand=self.region.strand,
                                name="intron"
                            )
                        )
                        start += l

                    elif c == 4:  # soft clip
                        continue

                    else:
                        continue

                features_list = []

                # for m6a
                m6a_loci = -1
                polya_length = -1
                if self.features:
                    if read.has_tag(self.features["m6a"]):
                        # multiple m6a support
                        m6a_loci_list = list(
                            map(float, [loci.strip() for loci in read.get_tag(self.features["m6a"]).split(',')])
                        )
                        for m6a_loci in m6a_loci_list:
                            features_list.append(GenomicLoci(
                                chromosome=self.region.chromosome,
                                start=m6a_loci,
                                end=m6a_loci,
                                strand="*",
                                name="m6a"
                            ))

                    if read.has_tag(self.features["real_strand"]) and read.has_tag(self.features["polya"]):
                        polya_length = float(read.get_tag(self.features["polya"]))
                        real_strand = read.get_tag(self.features["real_strand"])
                        assert read.get_tag(self.features["real_strand"]) in {"+", "-", "*"}, \
                            f"strand should be *, + or -, not {real_strand}"
                        if real_strand == "+":
                            features_list.append(
                                GenomicLoci(
                                    chromosome=self.region.chromosome,
                                    start=read.reference_end + 1,
                                    end=read.reference_end + polya_length - 1,
                                    strand=real_strand,
                                    name="polya"
                                )
                            )
                        elif real_strand == "-":
                            features_list.append(
                                GenomicLoci(
                                    chromosome=self.region.chromosome,
                                    start=read.reference_start - polya_length + 1,
                                    end=read.reference_start,
                                    strand=real_strand,
                                    name="polya"
                                )
                            )
                        else:
                            pass

                read = Reads(
                    chromosome=self.region.chromosome,
                    start=min(map(lambda x: x.start, exon_bound + features_list)),
                    end=max(map(lambda x: x.end, exon_bound + features_list)),
                    strand=self.region.strand,
                    id=read.query_name,
                    exons=exon_bound,
                    introns=intron_bound,
                    polya_length=polya_length,
                    m6a=m6a_loci,
                    features=features_list
                )

                if read.start < self.region.start or read.end > self.region.end:
                    continue

                self.data.append(read)

        except IOError as err:
            logger.error('There is no .bam file at {0}'.format(self.path))
            logger.error(err)
        except ValueError as err:
            logger.error(self.path)
            logger.error(err)

    def load(
            self,
            region: GenomicLoci,
            *args,
            **kwargs):
        u"""
        loading data
        :param region: the plotting region
        :return:
        """

        self.region = region

        # load from bed file
        if self.is_bed:
            self.load_bed()
        else:
            self.load_bam()

        tmp_df = pd.DataFrame(
            map(lambda x: x.to_dict(), self.data)
        )
        tmp_df["list_index"] = range(len(self.data))
        if self.exon_focus:
            e_f_use = []
            for read in self.data:
                tmp_ind = []
                for f_e in self.exon_focus:
                    tmp_ind.append("1" if f_e in str(read) else "0")
                e_f_use.append("_".join(tmp_ind))
            tmp_df["exon_group"] = e_f_use
        else:
            tmp_df["exon_group"] = "0"
        self.meta = self.df_sort(tmp_df)

    def len(self, scale: Union[int, float] = 0.005) -> int:
        u"""
        the length of reference to draw in final plots, default using the quarter of number of transcripts
        """
        size = len(self.data)
        return int(max(size * scale, 1))


if __name__ == '__main__':
    pass

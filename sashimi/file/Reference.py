#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2020.05.07

This scripts contains the class handle the reference file
"""
import glob
import gzip
import os
import re
from collections import namedtuple, defaultdict
from typing import List, Union, Optional

import filetype
import pysam
from loguru import logger

from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.Protein import CdsProtein
from sashimi.base.Readder import Reader
from sashimi.base.Transcript import Transcript
from sashimi.file.File import File


class Reference(File):
    u"""
    The reference file, support gtf and gff format
    """

    def __init__(self, path: str, category: str = "gtf",
                 add_domain: bool = False, add_local_domain: Optional[str] = False,
                 proxy: Optional[str] = None, timeout: int = 10):
        u"""
        init func
        :param path: path to input file
        """
        categories = ["gtf", "bam", "bed"]
        assert category in categories, f"category should be one of {categories}, instead of {category}"

        super().__init__(path=self.index_gtf(path) if category == "gtf" else path)
        self.category = category
        self.data = []
        self.add_domain = add_domain
        self.domain = None
        self.interval_file = {}

        self.add_local_domain = add_local_domain
        self.local_domain = None

        self.proxy = proxy
        self.timeout = timeout

        if proxy:
            logger.info(f"Using proxy: {proxy}")

    def __add__(self, other):
        assert isinstance(other, Reference), "only Reference and Reference could be added"
        new_ref = Reference(self.path, category=self.category)
        new_ref.data += self.data
        new_ref.data += other.data
        new_ref.data = sorted(new_ref.data)

        new_ref.interval_file.update(other.interval_file)

        if self.domain:
            new_ref.__add_domain__()
        if self.add_local_domain:
            new_ref.__load_local_domain__()
        return new_ref

    def len(self, scale: Union[int, float] = .25) -> int:
        u"""
        the length of reference to draw in final plots, default using the quarter of number of transcripts
        """
        size = len(self.data)

        if self.domain:
            size += len(self.domain)

        return int(max(size * float(scale), 1))

    @property
    def exons(self) -> List[List[int]]:
        res = []
        for transcript in self.data:
            for exon in transcript.exons:
                res.append([exon.start, exon.end])
        return sorted(res, key=lambda x: [x[0], x[1]])

    @classmethod
    def create(cls, path: str, add_domain: bool = False, add_local_domain: Optional[str] = False,
               category: str = "gtf"):
        u"""
        create reference file object
        :param path: path to input file
        :param add_domain: whether plot domain
        :param add_local_domain: fetch domain information from local file, which download from ucsc
        :param category: the type of reference file, include gtf and bam: customized reads as references
        :return: Reference obj
        """
        assert os.path.exists(path), f"{path} not exists"

        if add_local_domain:
            assert os.path.isdir(add_local_domain), f"{add_local_domain} not exists"

        return cls(
            path=path,
            add_domain=add_domain,
            add_local_domain=add_local_domain,
            category=category)

    def __load_local_domain__(self, region: GenomicLoci):

        if self.local_domain:
            u"""
            Because there is no gene or transcript id, we don't merge multiple domain here.
            """
            pass

        Domain_region = namedtuple(
            "DomainGenomicRegion",
            ["category", "type", "description", "unique_id", "start", "end"]
        )
        protein_info = defaultdict(list)

        if self.add_local_domain:
            for bb_file in glob.glob(pathname=f"{self.add_local_domain}/*.bb"):
                domain_res = defaultdict(list)
                base_name = os.path.basename(bb_file).replace(".bb", "")
                for record in Reader.read_bigbed(bb_file, region):
                    record = record[2].split("\t")
                    current_id = record[0]
                    current_start = int(record[3])
                    block_sizes = [int(x) for x in record[7].split(",") if x]
                    block_starts = [int(x) for x in record[8].split(",") if x]
                    current_desc = record[17]

                    current_domain_res = []
                    for index in range(len(block_starts)):
                        current_domain_res.append(
                            Domain_region._make(
                                [current_id,
                                 current_id,
                                 current_desc,
                                 current_id,
                                 current_start + 1 + block_starts[index],
                                 current_start + 1 + block_starts[index] + block_sizes[index] - 1
                                 ]
                            )
                        )
                    domain_res[current_id].extend([current_domain_res])

                for domain_unique_id, domain_list in domain_res.items():
                    start_site = min([
                        min(map(lambda x: x.start, i)) for i in domain_list
                    ])

                    end_site = max([
                        max(map(lambda x: x.end, i)) for i in domain_list
                    ])

                    protein_info[base_name].append(
                        Transcript(
                            chromosome=region.chromosome,
                            start=start_site,
                            end=end_site,
                            strand="*",
                            exons=domain_list,
                            gene=domain_list[0][0].unique_id,
                            domain_type=domain_list[0][0].type,
                            domain_description=domain_list[0][0].description,
                            domain_category=domain_list[0][0].category,
                            category="protein"
                        )
                    )

            self.local_domain = protein_info

    def __add_domain__(self):
        gene_id = set(map(lambda x: x.gene_id, self.data))
        transcript_id = set(map(lambda x: x.transcript_id, self.data))
        chromosome_id = self.data[0].chromosome

        # if domain is not None, then add it.
        if self.domain:
            gene_id = gene_id.difference(self.domain.gene_id)
            transcript_id = transcript_id.difference(self.domain.transcript_id)

        if len(transcript_id) != 0 and len(gene_id) != 0:
            with pysam.Tabixfile(self.path) as gtf_tabix:
                domain_info = CdsProtein.__re_iter_gtf__(
                    gtf_tabix=gtf_tabix,
                    chromosome=chromosome_id,
                    transcript_id=transcript_id,
                    gene_id=chromosome_id
                )
            # add pep information to cdsProtein.
            if self.domain:
                self.domain.add(domain_info)
            else:
                self.domain = domain_info

    def add_interval(self, interval: str, label: str):
        u"""
        Add another annotation in to data, and each annotation has a track
        :param interval: a bed file
        :param label
        :return:
        """
        assert os.path.exists(interval), f"{interval} not exists."
        self.interval_file[interval] = label
        return self

    @staticmethod
    def is_gtf(infile):
        u"""
        check if input file is gtf
        :param infile: path to input file
        :return:
        """
        if infile is None:
            return False

        is_gtf = 0
        try:
            if filetype.guess_mime(infile) == "application/gzip":
                is_gtf += 10
                r = gzip.open(infile, "rt")
            else:
                r = open(infile)

            for line in r:
                if line.startswith("#"):
                    continue

                lines = re.split(r"\s+", line)

                if len(lines) < 8:
                    break

                if re.search(r"([\w-]+ \"[\w.\s\-%,:]+\";? ?)+", " ".join(lines[8:])):
                    is_gtf += 1

                break

            r.close()
            return is_gtf
        except TypeError as err:
            logger.error("failed to open %s", infile)
            exit(err)

    @classmethod
    def sort_gtf(cls, input_gtf: str, output_gtf: str):

        if os.path.exists(output_gtf) and os.path.getsize(output_gtf) > 10:
            logger.info(f"{output_gtf} already exists, skip sorting")
            return

        data = []
        logger.info("Sorting %s" % input_gtf)

        try:
            w = open(output_gtf.replace(".gz", ""), "w+")
            r = gzip.open(input_gtf, "rt") if input_gtf.endswith(".gz") else open(input_gtf)
            for line in r:
                if line.startswith("#"):
                    w.write(line)
                    continue

                lines = line.split()

                if len(lines) < 1:
                    continue

                data.append(
                    GenomicLoci(
                        chromosome=lines[0], start=lines[3], end=lines[4],
                        strand=lines[6], gtf_line=line
                    )
                )
            r.close()

            for i in sorted(data, key=lambda x: [x.chromosome, x.start]):
                w.write(i.gtf_line)
            w.close()

            pysam.tabix_compress(output_gtf.replace(".gz", ""), output_gtf)
            os.remove(output_gtf.replace(".gz", ""))
        except IOError as err:
            logger.error(f"could not sort gtf, because {err}")
            exit(err)

    @classmethod
    def index_gtf(cls, input_gtf):
        u"""
        Check the tabix index of input gtf file

        Extract only exon tags and keep it clean
        :param input_gtf: path to input gtf file
        :return path to compressed and indexed bgzipped gtf file
        """
        gtf = cls.is_gtf(input_gtf)
        assert gtf % 10 == 1, f"{input_gtf} seems not be gtf format"

        if gtf // 10 > 0:
            output_gtf = input_gtf
        else:
            output_gtf = input_gtf + ".gz"

        if not os.path.exists(output_gtf + ".tbi"):
            logger.info(f"Create index for {input_gtf}")
            try:
                pysam.tabix_index(
                    output_gtf,
                    preset="gff",
                    force=True,
                    keep_original=True
                )
            except OSError as err:
                logger.error(err)
                logger.error("Guess gtf needs to be sorted")

                sorted_gtf = re.sub(r"\.gtf(.gz)?$", "", input_gtf) + ".sorted.gtf.gz"
                if os.path.exists(sorted_gtf) and os.path.exists(sorted_gtf + ".tbi"):
                    return sorted_gtf

                cls.sort_gtf(input_gtf, sorted_gtf)
                pysam.tabix_index(
                    sorted_gtf, preset="gff",
                    force=True, keep_original=True
                )
                return sorted_gtf
        elif os.path.getctime(output_gtf) < os.path.getctime(output_gtf):
            logger.info("the tbi index is older than the gtf file")

        return output_gtf

    def __load_gtf__(self, region: GenomicLoci) -> List[Transcript]:

        u"""
        Load transcripts inside of region from gtf file
        :param region: target region
        :return: list of Transcript
        """
        transcripts = {}
        exons = {}

        for rec in Reader.read_gtf(self.path, region):
            start = max(rec.start, region.start)
            end = min(rec.end, region.end)

            if end + 1 <= region.start:
                continue
            if start + 1 >= region.end:
                break

            if re.search(r"(rna|transcript|cds)", rec.feature, re.I):
                if rec.transcript_id not in transcripts.keys():
                    transcripts[rec.transcript_id] = Transcript(
                        chromosome=rec.contig,
                        start=start,
                        end=end,
                        strand=rec.strand,
                        transcript_id=rec.transcript_id,
                        gene_id=rec.gene_id,
                        gene=rec.gene_name if "gene_name" in rec.attributes else "",
                        transcript=rec.transcript_name if "transcript_name" in rec.attributes else "",
                        exons=[]
                    )
            elif re.search(r"(exon)", rec.feature, re.I):
                if rec.transcript_id not in exons.keys():
                    exons[rec.transcript_id] = []

                exons[rec.transcript_id].append(
                    GenomicLoci(
                        chromosome=rec.contig,
                        start=start,
                        end=end,
                        strand=rec.strand
                    )
                )

        for key, trans in transcripts.items():
            if key in exons.keys():
                trans.exons += exons[key]

        return sorted(transcripts.values())

    def __load_bam__(self, region: GenomicLoci, threshold_of_reads: int = 0) -> List[Transcript]:
        u"""
        Load transcripts inside of region from bam file
        :param region: target region
        :param threshold_of_reads: only kept reads with minimum frequency
        :return: list of Transcript
        """
        transcripts = {}
        try:
            for read, strand in Reader.read_bam(self.path, region):
                start = read.reference_start

                exons_in_read = []
                for cigar, length in read.cigartuples:
                    cur_start = start + 1
                    cur_end = start + length + 1

                    if cigar == 0:  # M
                        if cur_start < region.end < cur_end:
                            exons_in_read.append(GenomicLoci(
                                chromosome=read.reference_name,
                                start=cur_start if cur_start > region.start else region.start,
                                end=cur_end if cur_end <= region.end else region.end,
                                strand="+",
                            ))

                    elif cigar not in (1, 2, 4, 5):  # I, D, S, H
                        start += length

                t = Transcript(
                    chromosome=read.reference_name,
                    start=read.reference_start + 1 if read.reference_start + 1 > region.start else region.start,
                    end=read.reference_end + 1 if read.reference_end + 1 < region.end else region.end,
                    strand=strand,
                    exons=exons_in_read
                )
                transcripts[t] = transcripts.get(t, 0) + 1
        except IOError as err:
            logger.error('There is no .bam file at {0}'.format(self.path))
            logger.error(err)
        except ValueError as err:
            logger.error(self.path)
            logger.error(err)

        return sorted([x for x, y in transcripts.items() if y > threshold_of_reads])

    def __load_bed__(self, region: GenomicLoci) -> List[Transcript]:
        transcripts = []
        try:
            for rec in Reader.read_gtf(self.path, region=region, bed=True):
                exon_bound = []
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

                read = Transcript(
                    chromosome=self.region.chromosome,
                    start=min(map(lambda x: x.start, exon_bound)),
                    end=max(map(lambda x: x.end, exon_bound)),
                    strand=self.region.strand,
                    transcript_id=current_id,
                    exons=exon_bound,
                )
                if read.start < self.region.start or read.end > self.region.end:
                    continue

                transcripts.append(read)

        except IOError as err:
            logger.error('There is no .bed file at {0}'.format(self.path))
            logger.error(err)
        except ValueError as err:
            logger.error(self.path)
            logger.error(err)
        return transcripts

    def load(self, region: GenomicLoci, threshold_of_reads: int = 0, **kwargs):
        u"""
        Load transcripts inside of region
        :param region: target region
        :param threshold_of_reads: used for bam file, only kept reads with minimum frequency
        :return:
        """
        assert isinstance(region, GenomicLoci), "region should be a GenomicLoci object"
        self.region = region
        if self.category == "gtf":
            self.data = self.__load_gtf__(region)
            if self.add_local_domain:
                self.__load_local_domain__(region)

            if self.add_domain:
                self.__add_domain__()
        elif self.category == "bam":
            self.data = self.__load_bam__(region, threshold_of_reads)
        elif self.category == "bed":
            self.data = self.__load_bed__(region)

        rec, start, end, strand = None, None, None, None
        for interval_file, interval_label in self.interval_file.items():
            try:
                if not os.path.exists(interval_file + ".tbi"):
                    interval_file = pysam.tabix_index(
                        interval_file,
                        preset="bed",
                        force=True,
                        keep_original=True
                    )

                interval_target = []
                for rec in Reader.read_gtf(interval_file, region=self.region, bed=True):
                    start = max(rec.start, self.region.start)
                    end = min(rec.end, self.region.end)

                    if end + 1 <= self.region.start:
                        continue
                    if start + 1 >= self.region.end:
                        break

                    try:
                        strand = rec.strand
                        rec_name = rec.name
                    except KeyError:
                        strand = "*"
                        rec_name = ""

                    interval_target.append(
                        GenomicLoci(
                            chromosome=rec.contig,
                            start=rec.start,
                            end=rec.end,
                            strand=strand,
                            name=rec_name
                        )
                    )

                if len(interval_target) != 0 and rec is not None:
                    self.data.append(Transcript(
                        chromosome=rec.contig,
                        start=start,
                        end=end,
                        strand=strand,
                        exons=interval_target,
                        transcript_id=interval_label,
                        gene_id=interval_label,
                        gene=interval_label,
                        transcript=interval_label,
                        category="interval"
                    ))

            except NameError:
                raise "Target region was not found, run load before add_bed."

            except OSError:
                raise f"Error found when build index for {interval_file}, please sort file manually"


if __name__ == "__main__":
    pass

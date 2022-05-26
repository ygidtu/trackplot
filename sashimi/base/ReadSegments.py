#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Generate object for plotting reads like IGV track
"""
import os.path
from typing import Optional

import pandas as pd
import pysam

from conf.logger import logger
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.Readder import Reader
from sashimi.file.File import File


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
        :param features: current support m6a site and polya length, a list of GenomicLoci
        """

        super().__init__(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand
        )

        self.exons = sorted(exons)
        self.introns = sorted(introns)
        self.polya_length = polya_length
        self.m6a = m6a
        self.features = features
        self.id = id

    @property
    def exon_list(self):
        u"""
        return a nested list which contains exon regions
        :return: a nested list
        """

        exon_nested_lst = []
        for i in self.exons:
            exon_nested_lst.append(([i.start + 1, i.end]))
        return exon_nested_lst

    @property
    def intron_list(self):
        u"""
        return a nested list which contains intronic regions
        :return: a nested list
        """

        intron_nested_lst = []
        for i in self.introns:
            intron_nested_lst.append(
                ([i.start + 1, i.end])
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

    def to_dict(self):
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


class ReadSegment(File):

    def __init__(
            self,
            path: str,
            meta: Optional[pd.DataFrame] = None,
            region: Optional[GenomicLoci] = None,
            library: str = "fr-unstrand"
    ):
        u"""

        :param path:
        :param library:
        """
        super().__init__(path)

        assert library in ["fr-firststrand", "fr-secondstrand", "fr-unstrand"], \
            "Illegal library name."

        self.library = library
        self.data = []
        self.meta = meta
        self.region = region

    @classmethod
    def create(
            cls,
            path: str,
            library: str = "fr-unstrand"
    ):
        u"""

        :param path:
        :param library:
        :return:
        """
        if not os.path.exists(path + ".bai"):
            pysam.index(path)

        return cls(
            path=path,
            library=library
        )

    def set_region(self,
                   chromosome,
                   start,
                   end,
                   strand):
        self.region = GenomicLoci(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand
        )

    def load(
            self,
            region: GenomicLoci,
            features: Optional[dict] = None
    ):
        u"""
        Load each reads and its strand, features.
        m6a tag support the genomic site
        polya tag support length of polya. if polya tag is available, then read_strand (rl) is also essential
        :param region:
        :param features: support m6a and polyA length from bam tag. like {"m6a": "ma", "polya": "pa", "real_strand": "rs"}
        :return:
        """

        self.region = region
        feature_ids = ["m6a", "polya", "read_strand"]

        try:
            for read, _ in Reader.read_bam(self.path, self.region):

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
                                # strand information from minimap2
                                strand=self.region.strand,
                                name="exon"
                            )
                        )
                        start += l

                    elif c == 1:  # for insert
                        continue

                    elif c == 2:  # for del
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
                        # start += l

                    else:
                        continue

                features_list = []
                # for m6a
                m6a_loci = -1
                if read.has_tag(features["m6a"]):
                    m6a_loci = int(read.get_tag(features["m6a"]))
                    assert read.reference_start <= m6a_loci <= read.reference_end, \
                        f"{read.query_name}'s m6a loci was out of mapped region"

                    features_list.append(GenomicLoci(
                        chromosome=self.region.chromosome,
                        start=m6a_loci,
                        end=m6a_loci,
                        strand="*",
                        name="m6a"
                    ))

                polya_length = -1
                if features:
                    if read.has_tag(features["real_strand"]) and read.has_tag(features["polya"]):
                        polya_length = float(read.get_tag(features["polya"]))
                        real_strand = read.get_tag(features["real_strand"])
                        assert read.get_tag(features["real_strand"]) in {"+", "-", "*"}, \
                            f"strand should be *, + or -, not {real_strand}"
                        if real_strand == "+":
                            features_list.append(
                                GenomicLoci(
                                    chromosome=self.region.chromosome,
                                    start=read.reference_start + 1,
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
                                    end=read.reference_end,
                                    strand=real_strand,
                                    name="polya"
                                )
                            )
                        else:
                            pass

                self.data.append(
                    Reads(
                        chromosome=self.region.chromosome,
                        start=min(map(lambda x: x.start, exon_bound)),
                        end=min(map(lambda x: x.end, exon_bound)),
                        strand=self.region.strand,
                        id=read.query_name,
                        exons=exon_bound,
                        introns=intron_bound,
                        polya_length=polya_length,
                        m6a=m6a_loci,
                        features=features_list
                    )
                )

        except IOError as err:
            logger.error('There is no .bam file at {0}'.format(self.path))
            logger.error(err)
        except ValueError as err:
            logger.error(self.path)
            logger.error(err)

        self.meta = pd.DataFrame(
            map(lambda x: x.to_dict(), self.data)
        )

        self.meta["index"] = range(len(self.data))


if __name__ == '__main__':
    test = ReadSegment.create(path='../../example/bams/0.bam')

    test.load(
        region=GenomicLoci(
            chromosome="chr1",
            start=1270656,
            end=1284730,
            strand="+"),
        features={
            "m6a": "ma",
            "real_strand": "rs",
            "polya": "pa"
        })
    print(test.meta)

    test = ReadSegment.create(path='../../example/bams/WASH7P.bam')
    # chr1: 14362:29900
    test.load(
        region=GenomicLoci(
            chromosome="chr1",
            start=14362,
            end=29900,
            strand="+"),
        features={
            "m6a": "ma",
            "real_strand": "rs",
            "polya": "pa"
        })
    print(test.meta)

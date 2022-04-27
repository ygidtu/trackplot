#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""

from src.GenomicLoci import GenomicLoci


class Transcript(GenomicLoci):
    u"""
    Created by ygidtu at 2018.12.21

    A class inherit from GenomicLoci, to collect transcript information
    """

    __slots__ = [
        "gene",
        "gene_id",
        "transcript",
        "transcript_id",
        "exons",
        "is_reads",
        "show_id"
    ]

    def __init__(
        self, 
        chromosome: str, 
        start: int, 
        end: int,
        strand: str,
        exons: list, 
        gene: str = "",
        gene_id: str = "", 
        transcript: str = "",
        transcript_id: str = "", 
        is_reads: bool = False,
        show_id: bool = False
    ):
        u"""
        init this class
        :param chromosome: str
        :param start: int
        :param end: int
        :param strand: str
        :param gene_id: str
        :param exons: list of pysam.GTFProxy
        :param is_reads: is flag used by transcript  plot draw
        """

        super().__init__(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand
        )
        self.transcript = transcript
        self.transcript_id = transcript_id
        self.gene = gene
        self.gene_id = gene_id
        self.exons = exons
        self.is_reads = is_reads
        self.show_id = show_id

    def __str__(self):
        exons_str = []
        for i in self.exons:
            exons_str.append("{}-{}".format(i.start, i.end))

        return "{}:{}-{}:{} {} {}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand,
            self.transcript,
            "|".join(exons_str)
        )

    def __hash__(self):
        exons = sorted([str(x.__hash__()) for x in self.exons])
        return hash((self.chromosome, self.start, self.end, self.strand, " ".join(exons)))


if __name__ == "__main__":
    pass

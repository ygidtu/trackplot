#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06

Changelog:
    1. remove attributes
"""
from typing import List

from trackplot.base.GenomicLoci import GenomicLoci


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
        "category",
        "domain_category",
        "domain_type",
        "domain_description",
        "plot_intron"
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
            category: str = "exon",
            domain_category: str = "",
            domain_type: str = "",
            domain_description: str = ""
    ):
        u"""
        :param chromosome:
        :param start:
        :param end:
        :param strand:
        :param exons: A list of pysam.GTFProxy if category was exon, A nested tuple of list if category was protein
        :param gene: gene name when category is exon,such as "SAMD11"; domain's description when category is domain such as "Disordered"
        :param gene_id: gene id, such as "ENSG00000187634"
        :param transcript: transcript name, such as "SAMD11-011"
        :param transcript_id: transcript id, such as "ENST00000420190"
        :param category: exon or protein or interval
        :param domain_category: category of domain
        :param domain_description: description of domain
        :param domain_type: if category is protein, the type information of the given domain
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
        self.exons = sorted(exons)
        self.category = category
        self.domain_category = domain_category
        self.domain_type = domain_type
        self.domain_description = domain_description

    @property
    def exon_list(self):

        exon_nested_lst = []
        for i in self.exons:
            exon_nested_lst.append(([i.start + 1, i.end]))
        return exon_nested_lst

    def __str__(self):
        exons_str = []
        for i in self.exons:
            if isinstance(i, list):
                """
                2022.07.05
                Domain setting
                """
                exons_str.append("|".join(map(lambda x: f"{x.start}-{x.end}", i)))

            else:
                exons_str.append("{}-{}".format(i.start, i.end))

        return "{}:{}-{}:{} {} {}".format(
            self.chromosome,
            self.start,
            self.end,
            self.strand,
            self.transcript,
            "|".join(exons_str)
        )

    def __len__(self):
        return sum(map(lambda x: x[1] - x[0] + 1, self.exon_list))

    def __hash__(self):
        exons = sorted([str(x.__hash__()) for x in self.exons])
        return hash((self.chromosome, self.start, self.end, self.strand, " ".join(exons)))

    def __lt__(self, other):
        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome

        if self.start != other.start:
            return self.start < other.start

        if self.end != other.end:
            return self.end < other.end

        return len(self.exons) < len(other.exons)

    def __gt__(self, other):
        if self.chromosome != other.chromosome:
            return self.chromosome > other.chromosome

        if self.start != other.start:
            return self.start > other.start

        if self.end != other.end:
            return self.end > other.end

        return len(self.exons) > len(other.exons)

    def ids(self) -> List[str]:
        return [self.transcript, self.transcript_id, self.gene, self.gene_id]


if __name__ == "__main__":
    pass

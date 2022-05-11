#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""
from typing import List, Optional, Set
from copy import deepcopy

import numpy as np

from conf.logger import logger
from sashimi.anno.Stroke import Stroke
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.file.Reference import Reference
from sashimi.file.Fasta import Fasta
from sashimi.file.Bam import Bam
from sashimi.file.Bigwig import Bigwig


class Plot(object):
    u"""
    SpliceRegion represents a set of possible genomic segments
    this class is used to collect all the information about the exons and transcripts inside this region
    """

    def __init__(self):
        u"""
        init this class
        """
        self.region = None
        self.sites = {}
        self.focus = {}
        self.stroke = []
        self.events = None
        self.sequence = None
        self.reference = None
        self.graph_coords = None
        self.plots = []

    @property
    def chrom(self) -> Optional[str]:
        if self.region:
            return self.region.chrom

    @property
    def start(self) -> Optional[int]:
        if self.region:
            return self.region.start

    @property
    def end(self) -> Optional[int]:
        if self.region:
            return self.region.end

    def set_region(self, chromosome: str, start: int, end: int, strand: str = "+"):
        u"""
        change the plot region
        :param chromosome:
        :param start:
        :param end:
        :param strand:
        :return:
        """
        self.region = GenomicLoci(chromosome, start=start, end=end, strand=strand)

    def add_sites(self, sites):
        u"""
        highlight specific sites
        :param sites: string in 100,200 format or int
        :return:
        """
        if sites:
            if isinstance(sites, str):
                sites = [int(x) for x in sites.split(",")]

                for s in sites:
                    if s not in self.sites.keys():
                        self.sites[s - self.start] = "blue"
                    else:
                        self.sites[s - self.start] = "red"
            elif isinstance(sites, int):
                if sites not in self.sites.keys():
                    self.sites[sites - self.start] = "blue"
                else:
                    self.sites[sites - self.start] = "red"

    def add_focus(self, focus: Optional[str], start: int = 0, end: int = 0):
        u"""
        set focus region
        :param focus: string in 100-200:300-400
        :param start: start site
        :param end: end site
        :return:
        """
        if focus:
            for site in focus.split(":"):
                site = sorted([int(x) - self.start for x in site.split("-")])
                if site[0] < 0:
                    site[0] = 0

                if site[-1] > len(self):
                    site[-1] = len(self)

                self.focus[site[0]] = max(site[1], self.focus.get(site[0], -1))

        if 0 < start < end:
            self.focus[start] = max(end, self.focus.get(start, -1))

    def add_stroke(
            self,
            stroke: Optional[str],
            start: int = 0,
            end: int = 0,
            label: str = "",
            color: str = "black"
    ) -> List[Stroke]:
        u"""
        add stroke to plot
        :param stroke: stroke string in 100-200@red;300-400 format
        :param start: start site of stroke
        :param end: end site of stroke
        :param label: label of stroke
        :param color: color of stroke
        :return:
        """

        if stroke:
            self.stroke += Stroke.create_from_string(stroke)

        if 0 < start < end:
            self.stroke.append(Stroke(start, end, color, label))

    def set_sequence(self, fasta: str):
        u"""
        set sequence info for
        :param fasta: path to indexed fasta file
        :return:
        """
        logger.info(f"fetch sequence from {fasta}")
        self.sequence = Fasta.create(fasta)

    def set_reference(self, gtf: str):
        u"""
        add transcripts to this region
        :param gtf:
        :return:
        """
        logger.info(f"fetch transcripts from {gtf}")
        self.reference = Reference.create(gtf)

    def add_density(self,
                    path: str,
                    category: str = "bam",
                    label: str = "",
                    title: str = "",
                    barcodes: Optional[Set[str]] = None,
                    cell_barcode: str = "BC",
                    umi_barcode: str = "UB",
                    library: str = "fr-unstrand"
                    ):
        u"""
        add density object to plot
        :param label: the density y-axis label
        :param title: the density title
        :param barcodes: set of target barcodes, used by bam file of single cell sequencing
        :param cell_barcode: used by bam file of single cell sequencing
        :param umi_barcode: used by bam file of single cell sequencing
        :param library: the bam file library type
        :param path: the path to input file
        :param category: the input file type
        :return:
        """

        if category == "bam":
            obj = Bam.create(
                path,
                label=label,
                title=title,
                barcodes=barcodes,
                cell_barcode=cell_barcode,
                umi_barcode=umi_barcode,
                library=library
            )
        elif category == "bigwig" or category == "bw":
            obj = Bigwig.create(path, label=label)
        else:
            raise ValueError(f"the category should be one of [bam, bigwig, bw], instead of {category}")
        self.plots.append({"obj": obj, "type": "density"})

    def __len__(self) -> int:
        return self.end - self.start + 1

    def copy(self):
        return deepcopy(self)

    def init_graph_coords(self, read_depths_dict=None, intron_scale: float = 1, exon_scale: float = 1,
                          reverse_minus: bool = False):
        exon_starts = self.exon_starts
        exon_ends = self.exon_ends

        # modify information of reads
        if read_depths_dict:
            for key, val in read_depths_dict.items():
                for i in val.reads:
                    i.transcript = key.alias

                exon_starts += val.exon_starts
                exon_ends += val.exon_ends

        exon_coords = np.zeros(len(self))
        for i in range(len(exon_starts)):
            exon_coords[exon_starts[i] - self.start: exon_ends[i] - self.start] = 1

        graph_coords = np.zeros(len(self), dtype='f')

        x = 0
        for i in range(0, len(self)):
            graph_coords[i] = x

            if exon_coords[i if self.strand == '+' or not reverse_minus else -(i + 1)] == 1:
                x += 1. / exon_scale
            else:
                x += 1. / intron_scale

        self.graph_coords = graph_coords

    def get_relative(self, site: int) -> int:
        if self.graph_coords is None:
            raise ValueError("Please init graph_coords first")
        return self.graph_coords[site - self.start]


if __name__ == '__main__':
    pass

#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""
from typing import List, Optional
from copy import deepcopy

import numpy as np

from src.GenomicLoci import GenomicLoci
from src.Transcript import Transcript


class Stroke(object):
    def __init__(self, start: int, end: int, color: str = "red", label: str = ""):
        self.start = start
        self.end = end
        self.color = color
        self.label = label

    @property
    def center(self) -> float:
        return (self.end +  self.start) /  2


class SpliceRegion(GenomicLoci):
    u"""
    SpliceRegion represents a set of possible genomic segments
    this class is used to collect all the information about the exons and transcripts inside this region
    """

    def __init__(self, chromosome, start, end, strand, sites=None, events=None, ori=None, focus: Optional[str] = None,
                 stroke: Optional[str] = None):
        u"""
        init this class
        :param chromosome:  str, chromosome of these
        :param strand: str
        :param start: the very first site of list of exons
        :param end: the very last site of list of exons
        :param sites: list of int, all the splice sites
        :param events: the source splice events id
        :param focus: str -> 100-200:300-400, contains the focused regions
        """
        super().__init__(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand,
        )

        self.sites = self.set_sites(sites)
        self.focus = self.set_focus(focus)
        self.stroke = self.set_stroke(stroke)
        self.events = events
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.raster = False

        # {transcript_id: namedtuple(gtf proxy of transcript, [gtf proxy of exons])}
        self.__transcripts__ = {}
        self.ori = ori
        self.sequence = None

        self.__uniq_transcripts__ = set()

        # relative coords for plot
        self.graph_coords = None

    def set_sites(self, sites):
        res = {}

        if sites:
            for s in sites:
                if s not in res.keys():
                    res[s - self.start] = "blue"
                else:
                    res[s - self.start] = "red"
        return res

    def set_focus(self, focus: str):
        focus_sites = {}
        if focus:
            for site in focus.split(":"):
                site = sorted([int(x) - self.start for x in site.split("-")])
                if site[0] < 0:
                    site[0] = 0

                if site[-1] > len(self):
                    site[-1] = len(self)

                focus_sites[site[0]] = max(site[1], focus_sites.get(site[0], -1))
        return focus_sites

    def set_stroke(self, stroke: str) -> List[Stroke]:
        res = []

        if not stroke:
            return res

        for i in stroke.split(":"):
            i = i.split("@")
            sites = sorted([int(x) - self.start for x in i[0].split("-")])
            if sites[0] < 0:
                sites[0] = 0

            if sites[-1] > len(self):
                sites[-1] = len(self)

            color = "red"
            label = ""
            if len(i) > 1:
                i = i[-1].split("-")
                color = i[0]

                if len(i) > 1:
                    label = i[-1]

            res.append(Stroke(sites[0], sites[-1], color, label))

        return res

    def __str__(self):
        return '{0}:{1}-{2},{3}'.format(
            self.chromosome,
            self.start,
            self.end,
            self.__transcripts__
        )

    def __add__(self, other):
        u"""
        override add of genomic loci
        :param other: SpliceRegion
        :return:
        """
        if self.is_overlap(other):
            return SpliceRegion(
                chromosome=self.chromosome,
                start=min(self.start, other.start),
                end=max(self.end, other.end),
                strand=self.strand,
                sites=self.sites | other.sites
            )

    def __len__(self) -> int:
        return self.end - self.start + 1

    @property
    def exon_starts(self):
        u"""
        API for extract all exon starts
        :return:
        """
        starts = []
        for i in self.transcripts:
            for j in i.exons:
                starts.append(j.start)
        return sorted(starts)

    @property
    def exon_ends(self):
        u"""
        API for extract all exon ends
        :return:
        """
        ends = []
        for i in self.transcripts:
            for j in i.exons:
                ends.append(j.end)
        return sorted(ends)

    @property
    def exons(self):
        u"""
        API for extract all exons
        :return:
        """
        res = set()

        for i in self.transcripts:
            for j in i.exons:
                res.add(j)

        return sorted(res)

    @property
    def transcripts(self):
        u"""
        convert self.__transcripts__ to list format

        Note: wrapper gtf proxy to list and dict format
            1. there no need to change the code of sashimi plot
            2. shrink the transcript range

        :return: [[{transcript: id, gene: id, exon: []}, {}, {}], [{}]]
        """
        return sorted(
            [v for v in self.__transcripts__.values() if len(v.exons) > 0],
            key=lambda x: (x.start, x.end, len(x.exons)),
            reverse=True
        )

    def x_label(self, logtrans=None) -> str:
        label = f'Genomic coordinate ({self.chromosome}), "{self.strand}" strand'

        if logtrans is not None:
            label = f"{label}, y axis is log{logtrans} transformed"
        return label

    def add_gtf(self, gtf_line, show_id: bool = False):
        u"""
        add new gtf info to this Transcripts class
        :param gtf_line
        :param show_id: draw gene id and transcript id instead of gene name and transcript name
        :return:
        """
        if isinstance(gtf_line, Transcript):
            gtf_line = deepcopy(gtf_line)

            if gtf_line.start < self.start:
                gtf_line.start = self.start

            if gtf_line.end > self.end:
                gtf_line.end = self.end

            self.__transcripts__[gtf_line.transcript] = gtf_line

        else:
            if gtf_line.feature in ["transcript", "CDS"]:
                if gtf_line.transcript_id not in self.__transcripts__.keys():
                    self.__transcripts__[gtf_line.transcript_id] = Transcript(
                        chromosome=gtf_line.contig,
                        start=gtf_line.start + 1 if gtf_line.start + 1 > self.start else self.start,
                        end=gtf_line.end + 1 if gtf_line.end + 1 < self.end else self.end,
                        strand=gtf_line.strand,
                        transcript_id=gtf_line.transcript_id,
                        gene_id=gtf_line.gene_id,
                        gene=gtf_line.gene_name,
                        transcript=gtf_line.transcript_name,
                        exons=[],
                        show_id=show_id
                    )

            elif gtf_line.feature == "exon":
                if gtf_line.start + 1 >= self.end or gtf_line.end + 1 <= self.start:
                    return

                if gtf_line.transcript_id not in self.__transcripts__.keys():
                    raise ValueError("gtf file not sorted")

                tmp_transcript = self.__transcripts__[gtf_line.transcript_id]

                tmp_transcript.exons.append(
                    GenomicLoci(
                        chromosome=gtf_line.contig,
                        start=gtf_line.start + 1 if gtf_line.start + 1 > tmp_transcript.start else tmp_transcript.start,
                        end=gtf_line.end + 1 if gtf_line.end + 1 < tmp_transcript.end else tmp_transcript.end,
                        strand=gtf_line.strand
                    )
                )

    def get_region(self, genomic):
        u"""
        get smaller splice region
        :return:
        """

        tmp = SpliceRegion(
            chromosome=self.chromosome,
            start=genomic.start,
            end=genomic.end,
            strand=genomic.strand,
            events=genomic.events,
            sites=genomic.sites
        )

        for i in self.transcripts:
            if i.is_overlap(genomic):
                tmp.add_gtf(i)
        return tmp

    def filter_transcripts(self, transcripts: str):
        transcripts = transcripts.split(",")

        res = {}
        for k, v in self.__transcripts__.items():
            if v.transcript in transcripts or v.transcript_id in transcripts:
                res[k] = v
        self.__transcripts__ = res

    def copy(self):
        return deepcopy(self)

    def remove_empty_transcripts(self):
        u"""
        Remove transcript without exons
        :return:
        """
        res = {}
        for key, values in self.__transcripts__.items():
            if len(values.exons) > 0:
                res[key] = values

        self.__transcripts__ = res

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

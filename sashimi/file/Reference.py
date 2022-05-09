#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""
import gzip
import math
import os
import re

from typing import List, Optional

import filetype
import matplotlib as mpl
import numpy as np
import pysam

from matplotlib import pyplot as plt

from conf.logger import logger
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.Transcript import Transcript
from sashimi.anno.theme import Theme


class Reference(object):
    u"""
    The reference file, support gtf and gff format
    """

    def __init__(self, path: str):
        u"""
        init func
        :param path: path to input file
        """
        self.path = path
        self.transcripts = []

    def __add__(self, other):
        assert isinstance(other, Reference), "only Reference and Reference could be added"

        self.transcripts += other.transcripts
        self.transcripts = sorted(set(self.transcripts))

    @classmethod
    def create(cls, path: str):
        u"""
        create reference file object
        :param path: path to input file
        :return: Reference obj
        """
        assert os.path.exists(path), f"{path} not exists"
        return cls(path=cls.index_gtf(path))

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
        data = []
        logger.info("Sorting %s" % input_gtf)

        try:
            w = open(output_gtf.replace(".gz", ""), "w+")

            with open(input_gtf) as r:
                for line in r:
                    if line.startswith("#"):
                        fh.write(line)
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

        index = False
        if gtf // 10 > 0:
            output_gtf = input_gtf
        else:
            output_gtf = input_gtf + ".gz"

        sorted_gtf = re.sub(r"\.gtf(.gz)?$", "", input_gtf) + ".sorted.gtf.gz"
        if os.path.exists(sorted_gtf) and os.path.exists(sorted_gtf_gtf + ".tbi"):
            return sorted_gtf

        if not os.path.exists(output_gtf) or not os.path.exists(output_gtf + ".tbi"):
            logger.info(f"the {output_gtf} or tbi index not exists")
            index = True
        elif os.path.getctime(output_gtf) < os.path.getctime(output_gtf):
            logger.info("the tbi index is older than the gtf file")
            index = True

        if not index:
            return output_gtf

        logger.info("Create index for %s", input_gtf)
        try:
            cls.sort_gtf(input_gtf, sorted_gtf)
            pysam.tabix_index(
                input_gtf,
                preset="gff",
                force=True,
                keep_original=True
            )
        except OSError as err:
            if re.search("could not open", str(err)):
                raise err

            logger.error(err)
            logger.error("Guess gtf needs to be sorted")

            cls.sort_gtf(input_gtf, sorted_gtf)
            pysam.tabix_index(
                sorted_gtf, preset="gff",
                force=True, keep_original=True
            )

        return output_gtf

    def load(self, region: GenomicLoci):
        u"""
        Load transcripts inside of region
        :param region: target region
        :return:
        """
        assert isinstance(region, GenomicLoci), "region should be a GenomicLoci object"
        transcripts = {}

        r = pysam.TabixFile(self.path)

        try:
            iter_ = r.fetch(region.chromosome, region.start, region.end, parser=pysam.asGTF())
        except ValueError:
            try:
                if not region.chromosome.startswith("chr"):
                    logger.info("Guess need 'chr'")
                    iter_ = r.fetch(
                        "chr" + region.chromosome,
                        region.start, region.end, parser=pysam.asGTF()
                    )
                else:
                    logger.info("Guess 'chr' is redundant")
                    iter_ = r.fetch(
                        region.chromosome.replace("chr", ""),
                        region.start, region.end, parser=pysam.asGTF()
                    )
            except ValueError as err:
                logger.warn("please check the input region and gtf files")
                logger.error(err)
                raise err

        for rec in iter_:
            start = max(rec.start - region.start + 1, 0)
            end = max(rec.end - region.start + 1, 0)
            if re.search(r"(rna|transcript|cds)", rec.feature, re.I):
                if rec.transcript_id not in transcripts.keys():
                    transcripts[rec.transcript_id] = Transcript(
                        chromosome=rec.contig,
                        start=start,
                        end=end,
                        strand=rec.strand,
                        transcript_id=rec.transcript_id,
                        gene_id=rec.gene_id,
                        gene=rec.gene_name,
                        transcript=rec.transcript_name,
                        exons=[]
                    )

            elif re.search(r"(exon)", rec.feature, re.I):
                if rec.start + 1 >= region.end or rec.end + 1 <= region.start:
                    continue

                transcripts[rec.transcript_id].exons.append(
                    GenomicLoci(
                        chromosome=rec.contig,
                        start=start,
                        end=end,
                        strand=rec.strand
                    )
                )

        self.transcripts = sorted(transcripts.values())
        r.close()

    def plot(self,
             ax: mpl.axes.Axes,
             graph_coords: Optional[dict] = None,
             font_size: int = 5,
             distance_ratio: float = .3,
             show_gene: bool = False,
             show_id: bool = False,
             transcripts: Optional[List[str]] = None,
             remove_empty_transcripts: bool = False,
             color: Optional[str] = None,
             reverse_minus: bool = False,
             theme: str = "blank",
             y_loc: int = 0,
             exon_width: float = .3
             ):
        u"""
        draw the gene structure.

        :param ax: mpl.axes.Axes
        :param graph_coords: the convertor between genomic coords and plot coords
        :param font_size: the font size of transcript label
        :param show_gene: Boolean value to decide whether to show gene id in this plot
        :param show_id: Boolean value to decide whether to show the id or name of transcript/gene in this plot
        :param distance_ratio: distance between transcript label and transcript line
        :param transcripts: list of transcript id/name or gene id/name to specify which transcripts to show
        :param remove_empty_transcripts: do not show transcripts without any exons
        :param color: the color of transcripts
        :param reverse_minus:
        :param theme: plot style
        :param exon_width: scale of exons
        :param y_loc: default y-axis coords for labels
        :return:
        """
        logger.info("draw transcripts")
        Theme.set_theme(ax, theme)

        self.transcripts = sorted(self.transcripts)
        if not self.transcripts:
            return

        if not graph_coords:
            graph_coords = np.zeros(self.transcripts[-1].end - self.transcripts[0].start + 1, dtype=np.int)
            for i, j in enumerate(range(self.transcripts[0].start, self.transcripts[-1].end + 1)):
                graph_coords[j] = i

        """
        @2018.12.26
        Maybe I'm too stupid for this, using 30% of total length of x axis as the gap between text with axis
        """
        distance = distance_ratio * (max(graph_coords) - min(graph_coords))

        for transcript in self.transcripts:
            # ignore the unwanted transcript
            if transcripts and not (set(transcripts) & transcript.ids()):
                continue

            # ignore transcripts without any exons
            if remove_empty_transcripts and not transcript.exons:
                continue

            strand = transcript.strand
            # @2018.12.20 add transcript id, based on fixed coordinates
            if transcript.transcript:
                if show_gene and transcript.gene:
                    if show_id:
                        ax.text(
                            x=-1 * distance, y=y_loc + 0.25, s=transcript.gene_id,
                            fontsize=font_size
                        )

                        ax.text(
                            x=-1 * distance, y=y_loc - 0.25, s=transcript.transcript_id,
                            fontsize=font_size
                        )
                    else:
                        ax.text(
                            x=-1 * distance, y=y_loc, s=transcript.gene + " | " + transcript.transcript,
                            fontsize=font_size
                        )
                else:
                    ax.text(
                        x=-1 * distance, y=y_loc - 0.1, s=transcript.transcript,
                        fontsize=font_size
                    )

            # @2018.12.19
            # s and e is the start and end site of single exon
            for exon in transcript.exons:
                s, e, strand = exon.start, exon.end, exon.strand
                x = [
                    graph_coords[s], graph_coords[e],
                    graph_coords[e], graph_coords[s]
                ]
                y = [
                    y_loc - exon_width / 2, y_loc - exon_width / 2,
                    y_loc + exon_width / 2, y_loc + exon_width / 2
                ]
                ax.fill(x, y, 'k' if not color else color, lw=.5, zorder=20)

            # @2018.12.21
            # change the intron range
            # Draw intron.
            intron_sites = [graph_coords[transcript.start], graph_coords[transcript.end]]
            ax.plot(
                intron_sites, [y_loc, y_loc],
                color='k' if not color else color,
                lw=0.5
            )

            # @2018.12.23 fix intron arrows issues
            # Draw intron arrows.
            max_ = graph_coords[transcript.end]
            min_ = graph_coords[transcript.start]
            length = max_ - min_
            narrows = math.ceil(length / max(graph_coords) * 50)

            spread = .2 * length / narrows

            for i in range(narrows):
                loc = float(i) * length / narrows + graph_coords[transcript.start]
                if strand == '+' or reverse_minus:
                    x = [loc - spread, loc, loc - spread]
                else:
                    x = [loc + spread, loc, loc + spread]
                y = [y_loc - exon_width / 5, y_loc, y_loc + exon_width / 5]
                ax.plot(x, y, lw=.5, color='k')

            y_loc += 1  # if transcript.transcript else .5

        # ax.set_xlim(-2, max(graph_coords))
        ax.set_ylim(-.5, len(self.transcripts) + .5)


if __name__ == "__main__":
    pass

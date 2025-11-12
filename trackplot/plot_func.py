#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
This script contains the functions to draw different images
Modified by AD 2025/01/13
"""
import itertools # AD - for cumulative sum of graph_coords
#import gzip
import math
from copy import deepcopy
#from decimal import Decimal
from typing import Dict, List, Optional, Union, Tuple, Set

import matplotlib as mpl
import numpy as np
import seaborn as sns

from adjustText import adjust_text
from loguru import logger
from matplotlib import pylab
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch, Polygon
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import gaussian_kde, zscore
from scipy.stats import beta # AD

from trackplot.anno.theme import Theme
from trackplot.base.GenomicLoci import GenomicLoci
from trackplot.base.ReadDepth import ReadDepth
from trackplot.base.Stroke import Stroke
from trackplot.conf.config import CLUSTERING_METHOD, DISTANCE_METRIC
from trackplot.file.Annotation import Annotation
from trackplot.file.File import File
from trackplot.file.HiCMatrixTrack import HiCTrack
from trackplot.file.ReadSegments import ReadSegment


def load_barcodes(barcode: str) -> Tuple[Dict[str, Dict[str, Set[str]]], Dict[str, str]]:
    u"""
    as name says
    :param barcode: the path to barcode file
    """
    r = gzip.open(barcode, "rt") if barcode.endswith(".gz") else open(barcode, "r")
    res = {}
    colors = {}
    for line in r:
        line = line.strip().split()
        if len(line) >= 4:
            key, bc, group, color = line[:4]
            colors[group] = color
        elif len(line) == 3:
            key, bc, group = line[:3]
        elif len(line) == 2:
            key, bc = line
            group = ""
        else:
            continue

        if key not in res.keys():
            res[key] = {}

        if group not in res[key]:
            res[key][group] = set()

        res[key][group].add(bc)

    r.close()
    return res, colors


def get_limited_index(num, length):
    u"""
    Created by Zhang yiming at 2018.12.19
    Due to the original author didn't draw any element out of provided range
    So the scripts will through a lot of IndexError
    This function is used to scale that index into the reasonable range
    :param num: current index
    :param length: the list or numpy array length
    :return: (int, bool), 0 <= num <= length - 1, and modified or not
    """
    if num < 0:
        return 0, True

    if num >= length:
        return length - 1, True

    return num, False


def cubic_bezier(pts, t):
    """
    Get points in a cubic bezier.
    """
    p0, p1, p2, p3 = pts
    p0 = np.array(p0)
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    return p0 * (1 - t) ** 3 + 3 * t * p1 * (1 - t) ** 2 + 3 * t ** 2 * (1 - t) * p2 + t ** 3 * p3

def __merge_exons__(exons: List[List[int]]):
    u"""
    merge the overlap exons into one
    :param exons: sorted list of start and end sites of exons [[start, end], [start, end]]
    """
    exons = deepcopy(exons)
    res = []
    i, exon = 1, exons[0]
    while i < len(exons):
        curr_exon = exons[i]

        if curr_exon[0] <= exon[1]:
            exon[1] = curr_exon[1]
        else:
            res.append(exon)
            exon = curr_exon
        i += 1

    res.append(exon)
    return res


def init_graph_coords(region: GenomicLoci, exons: Optional[List[List[int]]] = None, exon_scale=1,
                      intron_scale=0.5) -> np.array:
    u"""
    init the default
    :param region: the plot region
    :param exons: list of start and end sites of exons [[start, end], [start, end]]
    :param exon_scale: the scale of exon, default set to 1
    :param intron_scale: the scale of intron, default set to 0.5 -> the intron showed in final will be half size of real
    """
    graph_coords = np.zeros(len(region), dtype=int)

    if exons:
        # AD - if user specifies an intron scale in line with Trackplot, addition is transparent
        if intron_scale <= 1:
            # if there have exons, init graph_coords by exon and intron scales
            for i in range(0, exons[0][0] - region.start):
                graph_coords[i] = (i - 0) * intron_scale
            exons = __merge_exons__(exons)
            for i in range(0, len(exons)):
                exon = exons[i]
                if i > 0:
                    intron = [exons[i - 1][1], exons[i][0]]

                    for j in range(intron[0], intron[1]):
                        if j >= region.start:
                            graph_coords[j - region.start] = graph_coords[intron[0] - region.start - 1] + (
                                    j - intron[0] + 1) * intron_scale
                for j in range(exon[0], exon[1] + 1):
                    if j >= region.start:
                        graph_coords[j - region.start] = graph_coords[exon[0] - region.start - 1] + (
                                j - exon[0] + 1) * exon_scale
            intron = [exons[-1][-1], region.end]
            for i in range(intron[0], intron[1]):
                if i >= region.start:
                    graph_coords[i - region.start] = graph_coords[intron[0] - region.start - 1] + (
                            i - intron[0] + 1) * intron_scale
        else:
            exons = __merge_exons__(exons)
            while exons[0][1] < region.start:
                exons = exons[1:]
            while exons[-1][0] > region.end:
                exons = exons[:-1]
            exons[0][0] = max(exons[0][0], region.start)
            exons[-1][1] = min(exons[-1][1], region.end)
            last_interval = region.end - exons[-1][1]
            for i, e in enumerate(exons):
                exons[i][0] -= region.start
                exons[i][1] -= region.start

            steps = [float(exon_scale)] * len(region)
            if exons[0][0] > 0:
                step = intron_scale / exons[0][0]
                for i in range(exons[0][0]):
                    steps[i] = step
            for e in range(1, len(exons)):
                interval = exons[e][0] - exons[e-1][1] - 2
                step = intron_scale / interval
                for i in range(exons[e-1][1], exons[e][0]):
                    steps[i] = step
                interval = exons[e][0] - exons[e-1][1] - 2
            if last_interval := len(region) - exons[-1][1] - 1:
                step = intron_scale / last_interval
                for i in range(exons[1][1] +1, len(region)):
                    steps[i] = step
            graph_coords = list(map(int, itertools.accumulate(steps)))
    else:
        # if there is not any exons, just init graph_coords by region
        for i, j in enumerate(range(region.start, region.end + 1)):
            graph_coords[j - region.start] = i

    # if there is no exon or intron in last of region
    if graph_coords[-1] == 0:
        current_max = max(np.where(graph_coords == np.max(graph_coords))[-1])
        for i in np.where(graph_coords == 0)[0]:
            if i > current_max:
                graph_coords[i] = max(graph_coords) + 1
    return graph_coords

def set_x_ticks(
        ax: mpl.axes.Axes,
        region: GenomicLoci,
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        sequence: Optional[Dict[int, str]] = None,
        log_trans: Optional[str] = None,
        nx_ticks: int = 4, font_size: int = 6, **kwargs):
    Theme.set_theme(ax, "blank")
    if graph_coords is None:
        graph_coords = init_graph_coords(region)

    # @2018.12.19 unnecessary text in figure
    x_label = str(region)

    if log_trans:
        if log_trans in ["2", "10"]:
            log_trans = f"log{log_trans}"
        x_label = f"{x_label}, y axis is {log_trans} transformed"

    ax.hlines(y=0, xmin=0, xmax=max(graph_coords), color="black", lw=1)
    ax.text(x=graph_coords[len(graph_coords) // 2], y=-2.8, s=x_label, fontsize=font_size, ha="center", va="top")

    # get x ticks coord by graph_coords
    bk = 1
    if not sequence and nx_ticks > 1:
        bk = max(graph_coords) // (nx_ticks - 1)

    # generate reverse coords map the genomic coords to x axis
    reverse_graph_coords = {}
    for i in range(1, len(graph_coords)):
        for x, y in zip(range(graph_coords[i-1], graph_coords[i]), np.linspace(i-1, i, num=graph_coords[i] - graph_coords[i-1])):
            reverse_graph_coords[x] = int(y)

    for i in range(graph_coords[0]):
        reverse_graph_coords[i] = graph_coords[0]

    line_space = {}
    for i in range(nx_ticks-1):
        i = bk * i
        line_space[i] = reverse_graph_coords[i] + region.start
    line_space[max(graph_coords)] = region.end

    if sequence:
        for i, seq in sequence.items():
            relative_i = graph_coords[i - region.start]
            if relative_i in line_space.keys():
                temp_txs = "{}\n{}".format(sequence.get(i, ""), line_space[relative_i])
            else:
                temp_txs = "\n{}".format(sequence.get(i, ""))
            line_space[relative_i] = temp_txs

    for x, s in line_space.items():
        ax.vlines(x=x, ymin=-.5, ymax=0, color="black", lw=1)
        ax.text(x=x, y=-1, s=s, fontsize=font_size, ha="center", va="top")

    ax.set_ylim(-3.5, 1)
    ax.set_xlim(min(graph_coords), max(graph_coords))


def set_y_ticks(
        ax: mpl.axes.Axes,
        label: str,
        graph_coords: Union[Dict, np.array],
        max_used_y_val: Union[int, float],
        min_used_y_val: Optional[Union[int, float]] = None,
        distance_between_label_axis: float = 0,
        n_y_ticks: int = 4,
        theme: str = "ticks",
        font_size: int = 5,
        show_y_label: bool = True,
        set_label_only: bool = False,
        y_axis_skip_zero: bool = True,
        **kwargs
):
    u"""
    The y ticks are formatted here
    @2019.03.31 add little check here to make sure the y-axis shows the real value
    """
    # set y ticks, y label and label
    Theme.set_theme(ax, theme)
    ax.set_xlim(min(graph_coords), max(graph_coords))
    if min_used_y_val is None:
        min_used_y_val, _ = ax.get_ylim()

    curr_y_tick_labels = []
    if not set_label_only:
        max_ = max_used_y_val
        plus = 0.2
        while max_ > 10:
            max_ /= 10
            plus /= 10

        plus = (max_used_y_val - min_used_y_val) * plus
        ax.set_ylim(min_used_y_val - plus, plus + max_used_y_val)
        ax.spines["left"].set_bounds(min_used_y_val, max_used_y_val)

        # assign number of ticks to minus and plus axis
        assign_ticks_y = [int(x / (abs(min_used_y_val) + max_used_y_val) * n_y_ticks) for x in
                          [abs(min_used_y_val), abs(max_used_y_val)]]
        universal_y_ticks = pylab.linspace(min_used_y_val, 0, assign_ticks_y[0] + 1)
        for i in pylab.linspace(0, max_used_y_val, assign_ticks_y[1] + 1):
            universal_y_ticks = np.append(universal_y_ticks, i)

        # add zero into strand-aware plot
        universal_y_ticks = np.unique(np.append(universal_y_ticks, 0))
        universal_y_ticks = sorted(universal_y_ticks)

        for lab in universal_y_ticks:
            curr_y_tick_labels.append(f"{int(lab)}")
            """
            if y_axis_skip_zero and lab == 0:
                # Exclude label for 0
                curr_y_tick_labels.append("")
            elif max_used_y_val < 1e-2:
                curr_y_tick_labels.append(f"{Decimal(lab):.2E}")
            else:
                curr_y_tick_labels.append(f"{lab:.1f}" if lab % 1 != 0 else f"{int(lab)}")
            """
        u"""
        @2019.01.04
        If there is no bam file, draw a blank y-axis 
        """
        ax.set_yticks(universal_y_ticks)
        ax.set_yticklabels(curr_y_tick_labels, fontsize=font_size)
        ax.yaxis.set_ticks_position('left')

    """
    Plot y labels
    @2018.12.20 using BAM label as y label
    @2019.01.04 change the standards of distance between y label and y-axis
    """

    if show_y_label:
        
        def __dynamic_distance__(distance_between_label_axis: float, label: str, scale: int = 100) -> float:
            if distance_between_label_axis != 0:
                return -distance_between_label_axis

            return max(0.01, math.ceil(len(label) / 10) * 10 / scale) * -1

        curr_y_tick_labels = sorted(curr_y_tick_labels, key=lambda x: len(x), reverse=True)
        ax.text(
            x=__dynamic_distance__(distance_between_label_axis, curr_y_tick_labels[0] if curr_y_tick_labels else "") * max(graph_coords),
            y=(max_used_y_val + min_used_y_val) / 2,
            s=label, fontsize=font_size, ha="right"
        )

    if max_used_y_val is not None and min_used_y_val is not None:
        ax.set_ylim(ymin=min_used_y_val, ymax=max_used_y_val)


def set_focus(
        ax: mpl.axes.Axes,
        graph_coords: Union[Dict, np.array],
        focus: Dict[int, int]
):
    for left, right in focus.items():
        try:
            left, right = graph_coords[left], graph_coords[right]
            fill_x = [left, right, right, left]

            y1, y2 = ax.get_ylim()
            fill_y = [y1, y1, y2, y2]
            ax.fill(fill_x, fill_y, alpha=0.1, color='grey')
        except IndexError as err:
            logger.debug("focus region is out of bound: " + str(err))


def set_indicator_lines(
        ax: mpl.axes.Axes,
        graph_coords: Union[Dict, np.array],
        sites: Dict[int, str],
        min_y_used: Union[int, float] = 0,
        max_y_used: Union[int, float] = None
):
    if sites is None:
        return

    if not max_y_used:
        min_y_used, max_y_used = ax.get_ylim()

    for site, color in sites.items():
        try:
            ax.vlines(
                x=graph_coords[site],
                ymin=min_y_used,
                ymax=max_y_used,
                color=color,
                linestyles="dashed",
                lw=0.5
            )
        except IndexError as err:
            logger.debug("Indicator line is out of bound: " + str(err))


def plot_stroke(
        ax: mpl.axes.Axes,
        data: List[Stroke],
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        font_size: int = 5,
        theme: str = "blank",
        **kwargs
):
    u"""
    plot stoke
    """

    strokes = sorted(data, key=lambda x: [x.start, x.end])

    for i, stroke in enumerate(strokes):
        try:
            ax.hlines(y=i, xmin=graph_coords[stroke.start], xmax=graph_coords[stroke.end],  color=stroke.color, lw=2)
            ax.text(-.1 * max(graph_coords),  i + .2, stroke.label, fontsize=font_size, color=stroke.color)
        except IndexError as _:
            logger.info(f"stroke {stroke} is out of bound, this stroke won't show up in final plot")

    Theme.set_theme(ax, theme)
    ax.set_xlim(min(graph_coords), max(graph_coords))
    ax.set_ylim(bottom=-1, top=len(strokes))


def plot_annotation(
        ax: mpl.axes.Axes,
        obj: Annotation,
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        font_size: int = 5,
        show_gene: bool = False,
        show_id: bool = False,
        color: Optional[str] = None,
        theme: str = "blank",
        y_loc: int = 0,
        exon_width: float = .3,
        plot_domain: bool = False,
        show_exon_id: bool = False,
        raster: bool = True,
        **kwargs
):
    u"""
    Plot the structure of annotation
    """
    Theme.set_theme(ax, theme)
    ax.set_xlim(min(graph_coords), max(graph_coords))

    region = obj.region

    data = sorted(obj.data)
    if not data:
        return

    if graph_coords is None:
        graph_coords = init_graph_coords(region)

    color = "k" if not color else color

    """
    @2018.12.26
    Maybe I'm too stupid for this, using 30% of total length of x axis as the gap between text with axis
    """
    patches, arrows = [], {}
    for transcript in data:
        strand = transcript.strand
        # @2018.12.20 add transcript id, based on fixed coordinates
        # @2025.02.04 add warnings for transcript name and id check
        if transcript.transcript or transcript.transcript_id:
            show_text = transcript.transcript
            if not transcript.transcript:
                logger.warning(f"there is not transcript_name, using transcript_id instead")
                show_text = transcript.transcript_id
            
            if show_gene and transcript.gene and transcript.gene_id != transcript.transcript_id:
                if show_id:
                    ax.text(x=-.01 * max(graph_coords), y=y_loc + 0.25,
                            s=transcript.gene_id, fontsize=font_size, ha="right")
                    ax.text(x=-.01 * max(graph_coords), y=y_loc - 0.25,
                            s=transcript.transcript_id, fontsize=font_size, ha="right")
                else:
                    ax.text(x=-.01 * max(graph_coords), y=y_loc,
                            s=transcript.gene + " | " + show_text, fontsize=font_size, ha="right")
            else:
                ax.text(x=-1, y=y_loc - 0.1, s=show_text, fontsize=font_size, ha="right")
        else:
            logger.warning(f"there is no transcript_name and transcript_id")

        # @2018.12.19
        # s and e is the start and end site of single exon
        # @2022.05.13
        # add index to avoid label overlapping of neighbor exon
        for ind, exon in enumerate(transcript.exons):
            s, e, strand = (region.relative(exon.start),
                            region.relative(exon.end), exon.strand)

            if s < 0 or e < 0:
                continue

            patches.append(plt.Rectangle(
                (graph_coords[s], y_loc - exon_width / 2),
                width=graph_coords[e]-graph_coords[s], height=exon_width,
                lw=.5, zorder=20,rasterized=raster, color=color
            ))

            if show_exon_id and exon.name:
                y_loc_offset = 0.1 if ind % 2 == 0 else - 0.2
                ax.text(x=(graph_coords[s] + graph_coords[e]) / 2,
                        y=y_loc + y_loc_offset,
                        s=exon.name, fontsize=font_size / 2,
                        ha="center", color="red")

        # @2018.12.21
        # change the intron range
        # Draw intron.
        # if transcript's category is interval, then don't plot the intron.
        if transcript.category != "interval":
            intron_sites = [
                graph_coords[region.relative(transcript.start)],
                graph_coords[region.relative(transcript.end)]
            ]
            ax.plot(intron_sites, [y_loc, y_loc], color=color, lw=0.5)

            # @2018.12.23 fix intron arrows issues
            # Draw intron arrows.
            max_ = graph_coords[region.relative(transcript.end)]
            min_ = graph_coords[region.relative(transcript.start)]
            length = max_ - min_
            narrows = math.ceil(length / max(graph_coords) * 50)

            if narrows > 0:
                spread = .2 * length / narrows

                for i in range(narrows):
                    loc = float(i) * length / narrows + graph_coords[region.relative(transcript.start)]

                    arrows[(loc - spread if strand == "+" else loc + spread, y_loc - exon_width / 5)] = Path.MOVETO
                    arrows[(loc, y_loc)] = Path.LINETO
                    arrows[(loc - spread if strand == "+" else loc + spread, y_loc + exon_width / 5)] = Path.LINETO

        y_loc += 1  # if transcript.transcript else .5
        if plot_domain and obj.domain and transcript.transcript_id in obj.domain.pep:
            current_domains = obj.domain.pep[transcript.transcript_id]

            for sub_current_domain in current_domains:
                if not (sub_current_domain.start <= region.end and
                        sub_current_domain.end >= region.start):
                    continue

                for sub_exon in sub_current_domain.exons:

                    for exon in sub_exon:
                        s, e = region.relative(
                            exon.start), region.relative(exon.end)
                        if e <= 0 or s > len(region):
                            continue
                        s = 0 if s < 0 else s
                        e = len(region) - 1 if e > len(region) else e

                        patches.append(plt.Rectangle(
                            (graph_coords[s], y_loc - exon_width / 2),
                            width=graph_coords[e] - graph_coords[s], height=exon_width,
                            lw=.5, zorder=20,rasterized=raster, color=color
                        ))

                    # @2022.05.13
                    intron_relative_s = region.relative(
                        min(map(lambda x_: x_.end, sub_exon)))
                    intron_relative_s = intron_relative_s if intron_relative_s >= 0 else 0
                    if intron_relative_s > len(region):
                        continue

                    intron_relative_e = region.relative(
                        max(map(lambda x_: x_.start, sub_exon)))
                    intron_relative_e = len(
                        region) - 1 if intron_relative_e > len(region) else intron_relative_e
                    if intron_relative_e <= 0:
                        continue

                    if len(sub_exon) != 1:
                        patches.append(plt.Rectangle(
                            (graph_coords[intron_relative_s], y_loc),
                            width=graph_coords[intron_relative_e] - graph_coords[intron_relative_s],
                            height=0, color=color, lw=0.2, rasterized=raster
                        ))
                if show_id:
                    ax.text(x=-1, y=y_loc - 0.125, s=f"{sub_current_domain.gene}|{transcript.transcript_id}",
                            fontsize=font_size / 2, ha="right")
                else:
                    ax.text(x=-1, y=y_loc - 0.125, s=f"{sub_current_domain.gene}|{transcript.transcript}",
                            fontsize=font_size / 2, ha="right")

                y_loc += .5

    if obj.local_domain:
        for base_name, current_domain in obj.local_domain.items():
            for sub_current_domain in current_domain:
                if not (sub_current_domain.start <= region.end and
                        sub_current_domain.end >= region.start):
                    continue

                for sub_exon in sub_current_domain.exons:

                    for exon in sub_exon:
                        s, e = region.relative(
                            exon.start), region.relative(exon.end)
                        if e <= 0 or s > len(region):
                            continue
                        s = 0 if s < 0 else s
                        e = len(region) - 1 if e > len(region) else e

                        patches.append(plt.Rectangle(
                            (graph_coords[s], y_loc - exon_width / 2),
                            width=graph_coords[e] - graph_coords[s], height=exon_width,
                            lw=.5, zorder=20, rasterized=raster, color=color
                        ))

                    # @2022.05.13
                    intron_relative_s = region.relative(
                        min(map(lambda x_: x_.end, sub_exon)))
                    intron_relative_s = intron_relative_s if intron_relative_s >= 0 else 0
                    if intron_relative_s > len(region):
                        continue

                    intron_relative_e = region.relative(
                        max(map(lambda x_: x_.start, sub_exon)))
                    intron_relative_e = len(
                        region) - 1 if intron_relative_e > len(region) else intron_relative_e
                    if intron_relative_e <= 0:
                        continue

                    if len(sub_exon) != 1:
                        patches.append(plt.Rectangle(
                            (graph_coords[intron_relative_s], y_loc),
                            width=graph_coords[intron_relative_e]-graph_coords[intron_relative_s],
                            height=0, color=color, lw=0.2, rasterized=raster
                        ))

                ax.text(x=-1, y=y_loc - 0.125, s=f"{sub_current_domain.gene}|{base_name}",
                        fontsize=font_size / 2, ha="right")

            y_loc += 1

    ax.add_collection(PatchCollection(patches, match_original=True))
    ax.add_patch(PathPatch(Path(list(arrows.keys()), list(arrows.values())), lw=.5, color=color, rasterized=raster))

    # @2022.05.13 Set y lim using y_loc value.
    ax.set_xlim(min(graph_coords), max(graph_coords))
    ax.set_ylim(-.5, y_loc + .5)


def plot_density(
        ax: mpl.axes.Axes,
        obj: Optional[File] = None,
        data: Optional[ReadDepth] = None,
        region: Optional[GenomicLoci] = None,
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        color="blue",
        font_size: int = 8,
        show_junction_number: bool = True,
        junction_number_font_size: int = 12,
        n_y_ticks: int = 4,
        distance_between_label_axis: float = .1,
        show_y_label: bool = True,
        y_label: str = "",
        theme: str = "ticks_blank",
        max_used_y_val: Optional[float] = None,
        min_used_y_val: Optional[float] = None,
        raster: bool = False,
        fill_step: str = "post",
        **kwargs
):
    u"""
    draw density plot
    :param ax: mpl.axes.Axes
    :param data: File
    :param data: ReadDepth
    :param region: GenomicLoci
    :param graph_coords:
    :param font_size: the font size for ticks, y-axis label and title
    :param show_junction_number: whether to show the number of junctions
    :param distance_between_label_axis: distance between y-axis label and y-axis ticks
    :param n_y_ticks: number of y ticks
    :param junction_number_font_size:
    :param obj: Bam or Bigwig object
    :param color: color for this density plot
    :param show_y_label: whether to show y-axis label
    :param y_label: the text of y-axis title
    :param theme: the theme name
    :param max_used_y_val: used to set same max y-axis
    :param raster:
    :param fill_step:
    :param kwargs:
    :return:
    """
    # max_used_y_val is None
    # distance_between_label_axis: 0
    if obj:
        assert obj.region is not None, "please load data first"
        region = obj.region
        if graph_coords is None:
            graph_coords = init_graph_coords(region)

        if not y_label:
            y_label = obj.label

    elif not (region and data):
        raise ValueError("please input obj or region and data")

    if data is None:
        data = obj.data

    try:
        jxns = data.junctions_dict
    except AttributeError:
        # depth do not have junctions
        jxns = {}
    
    # AD - convert junction counts to log scale early if requested
    if obj.log_trans and obj.log_trans.isdigit() and int(obj.log_trans) > 0: # in ["2", "10"]:

        y_label += f" (log{obj.log_trans})"
        denominator = np.log(int(obj.log_trans))
        for k, v in jxns.items():
            jxns[k] = np.log1p(v)  / denominator

    min_used_y_val = 0 # AD
    fixed_min_used_y = True #  AD
    fixed_max_used_y, fixed_min_used_y = max_used_y_val is not None, min_used_y_val is not None
    if max_used_y_val is None:
        if isinstance(data, dict):
            for k, v in data.items():
                max_used_y_val = max(v.plus.max(), max_used_y_val if max_used_y_val else 0)
        else:
            max_used_y_val = max(data.plus)
        if max_used_y_val % 2 == 1:
            max_used_y_val += 1

    if min_used_y_val is None:
        if isinstance(data, dict):
            min_used_y_val = [v.minus.max() if v.minus else 0 for v in data.values()]
            min_used_y_val = -1 * max(min_used_y_val)
        else:
            min_used_y_val = -1 * max(data.minus) if data.minus is not None else 0

    # Reduce memory footprint by using incremented graph_coords.
    x, y1, y2 = [], [], []
    for i in range(len(graph_coords)):
        x.append(graph_coords[i])

        if isinstance(data, dict):
            y1.append(max([v.plus[i] if v.plus is not None else 0 for v in data.values()]))
            y2.append(min([-v.minus[i] if v.minus is not None else 0 for v in data.values()]))
        else:
            y1.append(data.plus[i] if data.plus is not None else 0)
            y2.append(-data.minus[i] if data.minus is not None else 0)

    ax.fill_between(x, y1, y2=y2, color=color, lw=0, step=fill_step, rasterized=raster)

    if not isinstance(data, dict):
        if data.strand_aware:
            max_used_y_val = max(abs(min_used_y_val), max_used_y_val)
            min_used_y_val = -max(abs(min_used_y_val), max_used_y_val) if data.minus is not None else 0
            
    if jxns:
        # sort the junctions by intron length for better plotting look
        jxns_sorted_list = sorted(jxns.keys(), key=lambda x: (x.end - x.start, x.start, x.end), reverse=True)

        if not jxns:
            max_junction_count, min_junction_count = 0, 0
        else:
            max_junction_count = max(jxns.values())
            min_junction_count = min(jxns.values())
        junction_count_gap = max_junction_count - min_junction_count

        jxn_numbers = []
        for jxn_idx, jxn in enumerate(jxns_sorted_list):
            leftss, rightss = jxn.start, jxn.end

            # junction must at least have one anchor located in plotted region
            if (leftss < region.start < region.end < rightss) or rightss <= region.start or leftss >= region.end:
                logger.debug(f"junction {jxn} of {y_label} is is out of plotting region, skip")
                continue

            # @2018.12.19
            # set junctions coordinate here
            # the junction out of boundaries, set the boundaries as coordinate
            ss1_idx, ss1_modified = get_limited_index(leftss - region.start, len(graph_coords))
            ss2_idx, ss2_modified = get_limited_index(rightss - region.start, len(graph_coords))
            u"""
            @2019.01.14
            add two new variables to make it clear which one is index, which one is genomic site 
            """
            ss1, ss2 = graph_coords[ss1_idx], graph_coords[ss2_idx]
            # AD = keep junction arcs on top
            jxn_on_top = True

            # draw junction on bottom 
            if kwargs.get("density_by_strand"):
                jxn_on_top = jxn.strand == "+"
            else:
                jxn_on_top = jxn_idx % 2 == 0
                #jxn_on_top = True # AD - keep all junctions on same strand
                if abs(min_used_y_val) < max_used_y_val:
                    min_used_y_val = -max_used_y_val
            
            if fill_step == "post":
                ss1_idx = max(0, ss1_idx - 1)
            elif fill_step == "pre":
                ss2_idx = min(ss2_idx + 1, len(graph_coords) - 1)

            if jxn_on_top:
                left_dens, right_dens = data.curr_max(ss1_idx), data.curr_max(ss2_idx)
                current_height = abs(3 * max_used_y_val / 4)
                pts = [
                    left_dens if not ss1_modified else left_dens + current_height,
                    left_dens + current_height,
                    right_dens + current_height,
                    right_dens if not ss2_modified else right_dens + current_height
                ]
            else:
                left_dens, right_dens = abs(data.curr_min(ss1_idx)), abs(data.curr_min(ss2_idx))
                current_height = abs(3 * min_used_y_val / 4)
                pts = [
                    -left_dens if not ss1_modified else -left_dens - current_height,
                    -left_dens - current_height,
                    -right_dens - current_height,
                    -right_dens if not ss2_modified else -right_dens - current_height
                ]
  
            """
            @2018.12.26
            scale the junctions line width
            """
            if junction_count_gap > 0:
                #line_width = (jxns[jxn] - min_junction_count) / junction_count_gap
                line_width = .5 + 1.5 * (jxns[jxn] - min_junction_count) / junction_count_gap
            else:
                line_width = 0
            #line_width = max(.5,np.log())

            pts = [(ss1, pts[0]), (ss1, pts[1]), (ss2, pts[2]), (ss2, pts[3])]
            ax.add_patch(PathPatch(Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                                   ec=color, lw=line_width + 0.2, fc='none'))
            
            if show_junction_number:
                midpt = cubic_bezier(pts, .5)

                t = ax.text(
                    midpt[0], midpt[1],
                    '{0}'.format(round(jxns[jxn], 2)),
                    fontsize=junction_number_font_size,
                    ha='center', va='center',
                    backgroundcolor='w'
                )

                # @2018.12.19 transparent background
                t.set_bbox(dict(alpha=0))
                jxn_numbers.append(t)


    if obj and obj.title:
        ax.text(max(graph_coords) - len(obj.title), max_used_y_val, obj.title, color=color, fontsize=font_size)

    # update the y-axis actually used
    if not fixed_max_used_y:
        _, max_used_y_val = ax.get_ylim()

    if not fixed_min_used_y:
        min_used_y_val, _ = ax.get_ylim()

    if max_used_y_val == min_used_y_val == 0:
        min_used_y_val, max_used_y_val = ax.get_ylim()

    if not isinstance(data, dict) and data.strand_aware and kwargs.get("density_by_strand"):
        max_used_y_val = max(abs(min_used_y_val), max_used_y_val)
        min_used_y_val = -max_used_y_val
    elif not kwargs.get("density_by_strand") and not jxns:
        # if there is no any junctions, then just use the compressed y-axis
        min_used_y_val = 0

    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=max_used_y_val, # * 1.1,  # keep a buffer at top and bottom of junctions
        min_used_y_val=min_used_y_val, # * 1.1, # AD
        n_y_ticks=n_y_ticks,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        y_axis_skip_zero=False #if not isinstance(data, dict) and data.strand_aware else True
    )


def plot_site_plot(
        ax: mpl.axes.Axes,
        obj: File,
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        color="blue",
        font_size: int = 8,
        n_y_ticks: int = 3,
        distance_between_label_axis: float = .1,
        show_y_label: bool = True,
        y_label: str = "",
        strand_choice: str = None,
        theme="ticks",
        raster: bool = False,
        **kwargs
):
    """
    :param ax:
    :param obj:
    :param graph_coords:
    :param color:
    :param font_size:
    :param n_y_ticks:
    :param distance_between_label_axis:
    :param show_y_label:
    :param y_label:
    :param strand_choice:
    :param theme:
    :param raster:
    :return:
    """
    region = obj.region
    if graph_coords is None:
        graph_coords = init_graph_coords(region)

    if not y_label:
        y_label = obj.label

    if not show_y_label:
        y_label = ""

    plus, minus = obj.data.site_plus, obj.data.site_minus
    max_height, min_height = max(plus), min(minus)
    max_val = max(max_height, abs(min_height))

    patches = []
    for label, array_plot in zip(['plus', 'minus'], [plus, minus]):
        if strand_choice != "all" and label != strand_choice:
            continue

        array_hist = np.repeat(graph_coords, np.abs(array_plot).astype(np.int32))
        try:
            kde = gaussian_kde(array_hist)
            fit_value = kde.pdf(graph_coords)
        except (ValueError, np.linalg.LinAlgError):
            # logger.debug(err)
            # logger.debug(traceback.format_exc())
            continue

        fit_value = fit_value / fit_value.max()

        for x, y in zip(graph_coords, array_plot):
            patches.append(plt.Rectangle((x, 0), height=y, width=.6, color=color, rasterized=raster))
        ax.plot(graph_coords,
                fit_value * array_plot.max() if label == 'plus' else fit_value * array_plot.min(),
                c=color, lw=1)

    ax.add_collection(PatchCollection(patches, match_original=True))
    # set the y limit
    # set y ticks, y label and label
    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=1.1 * max_val,
        min_used_y_val=-1.1 * max_val,
        n_y_ticks=n_y_ticks,
        font_size=font_size,
        show_y_label=False,
        distance_between_label_axis=distance_between_label_axis,
        y_axis_skip_zero=False
    )


def plot_heatmap(
        ax: mpl.axes.Axes,
        cbar_ax: mpl.axes.Axes,
        data: Dict[str, ReadDepth],
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        color="viridis",
        font_size: int = 8,
        distance_between_label_axis: float = .1,
        show_y_label: bool = True,
        theme: str = "ticks",
        y_label: str = "",
        do_scale: bool = False,
        clustering: bool = False,
        clustering_method: str = "ward",
        distance_metric: str = "euclidean",
        raster: bool = True,
        show_row_names: bool = False,
        vmin=None, vmax=None,
        **kwargs
):
    u"""

    :param theme:
    :param graph_coords:
    :param ax: whether to scale the matrix
    :param cbar_ax: whether to scale the matrix
    :param font_size:
    :param data:
    :param distance_between_label_axis:
    :param show_y_label:
    :param y_label:
    :param do_scale: whether to scale the matrix
    :param clustering: whether reorder matrix by clustering
    :param clustering_method: same as  scipy.cluster.hierarchy.linkage
    :param distance_metric: same as scipy.spatial.distance.pdist
    :param color: used for seaborn.heatmap, see: https://matplotlib.org/3.5.1/tutorials/colors/colormaps.html
                'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
                'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'
    :param raster: whether to draw image in raster mode
    :param show_row_names:
    :param vmin: Values to anchor the colormap, otherwise they are inferred from the data and other keyword arguments.
    :param vmax: Values to anchor the colormap, otherwise they are inferred from the data and other keyword arguments.
    """
    labels = list(data.keys())
    mtx = np.array([x.wiggle for x in data.values()])

    if clustering and len(mtx) > 1:
        assert clustering_method in CLUSTERING_METHOD, f"clustering_method {clustering_method} is not supported."
        assert distance_metric in DISTANCE_METRIC, f"distance_metric {distance_metric} is not supported"

        order = dendrogram(
            linkage(mtx, method=clustering_method, metric=distance_metric),
            orientation='right',
            no_plot=True
        )

        mtx = mtx[order["leaves"], :]
        labels = [labels[x] for x in order["leaves"]]

    if do_scale:
        """
        y = (x – mean) / standard_deviation
        """
        mtx = zscore(mtx, axis=1)

    if not show_row_names:
        labels = False

    if vmin == vmax == 0 and (mtx.max() != 0 or mtx.min() != 0):
        logger.debug("The vmin and vmax == 0, but input matrix not != 0")
        vmin, vmax = None, None

    sns.heatmap(mtx, ax=ax, cmap=color, cbar_ax=cbar_ax,
                xticklabels=False, yticklabels=labels,
                center=False, rasterized=raster, vmin=vmin, vmax=vmax)

    ax.tick_params(axis='both', which='major', labelsize=font_size, rotation=0)
    cbar_ax.tick_params(labelsize=font_size)

    ymax, ymin = ax.get_ylim()
    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=ymax, min_used_y_val=ymin,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        set_label_only=True
    )


def plot_hic(
        ax: mpl.axes.Axes,
        cbar_ax: mpl.axes.Axes,
        obj: List[HiCTrack],
        show_legend: bool = True,
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        color: str = "RdYlBu_r",
        font_size: int = 8,
        distance_between_label_axis: float = .1,
        show_y_label: bool = True,
        theme: str = "ticks",
        y_label: str = "",
        n_y_ticks: int = 4,
        raster: bool = True,
        **kwargs
):
    assert len(obj) == 1, "HiC plot only support one file"

    obj = obj[0]

    if graph_coords is None:
        graph_coords = init_graph_coords(obj.region)

    if not y_label:
        y_label = obj.label

    y_max = obj.matrix.shape[1]
    ax.pcolormesh(obj.x - obj.region.start, obj.y, np.flipud(obj.matrix),
                  cmap=color, rasterized=raster)

    if obj.tad:

        for tad in obj.tad_list:
            x1 = tad.start - obj.region.start
            x2 = x1 + float(tad.end - tad.start) / 2
            x3 = tad.end - obj.region.start
            y1 = 0
            y2 = (x3 - x1)
            if y2 > obj.depth:
                continue

            triangle = Polygon([[x1, y1], [x2, y2], [x3, y1]], closed=False,
                               edgecolor="grey",
                               facecolor='None')

            ax.add_artist(triangle)

    ax.set_xlim(0, len(obj.region))
    ax.set_ylim(0, y_max)

    if show_legend:
        cbar = pylab.colorbar(
            mappable=mpl.cm.ScalarMappable(
                norm=mpl.colors.Normalize(
                    vmin=0, vmax=np.percentile(obj.matrix.diagonal(1), 80)),
                cmap=color),
            cax=cbar_ax)
        cbar.ax.tick_params(labelsize=font_size)
        if obj.log_trans:
            legend_ticks = cbar.get_ticks().tolist()
            legend_ticks[0] = f"{legend_ticks[0]}\n{obj.log_trans}"
            cbar.set_ticklabels(legend_ticks)

    ax.axis("off")
    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=obj.depth,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        n_y_ticks=n_y_ticks
    )


def plot_line(
        ax: mpl.axes.Axes,
        data: Dict[str, ReadDepth],
        graph_coords: Union[Dict, np.ndarray],
        font_size: int = 8,
        distance_between_label_axis: float = .1,
        show_y_label: bool = True,
        y_label: str = "",
        line_attrs: Optional[Dict[str, Dict]] = None,
        theme: str = "ticks_blank",
        n_y_ticks: int = 4,
        show_legend: bool = False,
        legend_position: str = "upper right",
        legend_ncol: int = 0,
        max_used_y_val: Optional[float] = None,
        **kwargs
):
    u"""

    :param ax: whether to scale the matrix
    :param font_size:
    :param data:
    :param distance_between_label_axis:
    :param show_y_label:
    :param graph_coords:
    :param y_label:
    :param line_attrs: dict contains attributes like color and type for each line
    :param theme: the theme name
    :param n_y_ticks: the number of y ticks
    :param show_legend: whether to show legend
    :param legend_position:
    :param legend_ncol:
    :param max_used_y_val:
    """
    max_y_val = 0
    legend = []
    distance_to_zero = pylab.linspace(0, len(graph_coords), len(data) + 2)
    for idx, (ylab, val) in enumerate(data.items()):
        attr = line_attrs.get(ylab, {}) if line_attrs else {}

        x, y = [], []
        for i, j in enumerate(val.wiggle):
            x.append(graph_coords[i])
            y.append(j)
        ax.plot(x, y, label=ylab if show_y_label else "", **attr)

        max_y_val = max(max_y_val, max(val.wiggle))

        if not show_legend:
            idx = int(distance_to_zero[idx + 1])
            t = ax.text(
                x[idx], np.median(y), ylab, color=attr.get("color", "black"),
                fontsize=font_size, ha='center', va='center', backgroundcolor='w'
            )

            # @2018.12.19 transparent background
            t.set_bbox(dict(alpha=0))
            legend.append(t)

    if show_legend:
        ax.legend(loc=legend_position,
                  ncol=int(len(data) / 1.5) if legend_ncol <= 0 else legend_ncol,
                  fancybox=False, shadow=False)
    else:
        try:
            adjust_text(legend, arrowprops=dict(arrowstyle="-", lw=1))
        except IndexError as err:
            logger.debug(err)

    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=max_y_val if max_used_y_val is None else max_used_y_val,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        n_y_ticks=n_y_ticks, **kwargs
    )


def plot_igv_like(
        ax: mpl.axes.Axes,
        obj: Dict[str, ReadSegment],
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        y_label: str = "",
        exon_color: Optional[str] = None,
        intron_color: Optional[str] = None,
        feature_color: Optional[str] = None,
        exon_width: float = .3,
        font_size: int = 8,
        n_y_ticks: int = 1,
        distance_between_label_axis: float = .1,
        show_y_label: bool = True,
        theme: str = "ticks_blank",
        raster: bool = False,
        **kwargs
):
    assert len(obj) == 1, "IGV-like plot only support one file"

    obj = list(obj.values())[0]
    assert isinstance(obj, ReadSegment), "the input obj should be ReadSegment"
    assert obj.region is not None, "please load data first"

    region = obj.region
    if graph_coords is None:
        graph_coords = init_graph_coords(region)

    if not y_label:
        y_label = obj.label
    y_loc = 0.5

    patches, scatters_x, scatters_y = [], [], []
    # Add this to skip zero coverage regions.
    if obj.meta is not None:

        for c_ind_list in obj.get_index():
            add_plot = False
            for c_ind in c_ind_list:
                c_data = obj.data[c_ind]
                # skip truncated reads in the given region
                if c_data.start < region.start or c_data.end > region.end:
                    continue

                for exon in c_data.exons:
                    s, e, strand = region.relative(
                        exon.start), region.relative(exon.end), exon.strand
                    if e < 0 or s > len(region):
                        continue

                    s = 0 if s < 0 else s
                    e = len(region) - 1 if e >= len(region) else e

                    patches.append(plt.Rectangle((graph_coords[s], y_loc - exon_width),
                                             width=graph_coords[e] - graph_coords[s],
                                             height=exon_width * 2,
                                             facecolor='k' if not exon_color else exon_color,
                                             lw=.5, zorder=20))

                    add_plot = add_plot | True

                for intron in c_data.introns:
                    s, e, strand = region.relative(
                        intron.start), region.relative(intron.end), intron.strand

                    if e < 0 or s > len(region):
                        continue
                    s = 0 if s < 0 else s
                    e = len(region) - 1 if e >= len(region) else e
                    intron_sites = [
                        graph_coords[s],
                        graph_coords[e]
                    ]
                    patches.append(plt.Rectangle(
                        (graph_coords[s], y_loc),
                        width=graph_coords[e]-graph_coords[s], height=0,
                        color="#4d4d4d" if not intron_color else intron_color,
                        lw=.2, rasterized=raster
                    ))

                for feature in c_data.features:
                    if feature.start == -1:
                        continue

                    s, e, strand = region.relative(feature.start), region.relative(
                        feature.end), feature.strand
                    is_site = False

                    if s == e:
                        s = s - 1
                        is_site = True

                    if e < 0 or s > len(region):
                        continue

                    s = 0 if s < 0 else s
                    e = len(region) - 1 if e >= len(region) else e

                    width_ratio = 0.75 if not is_site else 1.5

                    if is_site:
                        patches.append(plt.Rectangle((graph_coords[s], y_loc - exon_width * width_ratio),
                                                 width=graph_coords[e] - graph_coords[s],
                                                 height=exon_width * 2,
                                                 facecolor='blue' if not feature_color else feature_color,
                                                 lw=.2, zorder=20, rasterized = raster))

                        scatters_x.append(graph_coords[s])
                        scatters_y.append(y_loc + exon_width * width_ratio)
                    else:
                        patches.append(plt.Rectangle((graph_coords[s], y_loc - exon_width * width_ratio),
                                                 width=graph_coords[e] - graph_coords[s],
                                                 height=exon_width * 2,
                                                 facecolor='red' if not feature_color else feature_color,
                                                 lw=.2, zorder=20, rasterized = raster))
            if add_plot:
                y_loc += 1

    # using patches and draw scatter in bath to reduce time
    ax.add_collection(PatchCollection(patches, match_original=True))
    ax.scatter(
        scatters_x, scatters_y, rasterized=raster,
        color='b' if not feature_color else feature_color, s=.5, linewidths=(0,)
    )

    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=y_loc,
        n_y_ticks=n_y_ticks,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
    )


def plot_links(ax: mpl.axes.Axes,
               data: List[Stroke],
               graph_coords: Optional[Union[Dict, np.ndarray]] = None,
               max_y: int = -10, **kwargs):
    for stroke in sorted(data):
        leftss, rightss = graph_coords[stroke.start], graph_coords[stroke.end]

        step = (rightss - leftss) / 4

        pts = [
            (leftss, 0),
            (leftss + step, max_y),
            (rightss - step, max_y),
            (rightss, 0)
        ]
        a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        ax.add_patch(PathPatch(a, fc="none", ec=stroke.color, lw=1))

    ax.set_xlim(min(graph_coords), max(graph_coords))
    ax.set_ylim(max_y, 0)
    ax.axis("off")


def make_text_elements(text, x=0.0, y=0.0, width=1.0, height=1.0,
                       color='blue', edgecolor="black",
                       font=FontProperties(family='monospace')):
    tp = TextPath((0.0, 0.0), text, size=1, prop=font)
    bbox = tp.get_extents()
    bwidth = bbox.x1 - bbox.x0
    bheight = bbox.y1 - bbox.y0
    trafo = Affine2D()
    trafo.translate(-bbox.x0, -bbox.y0)
    trafo.scale(1 / bwidth * width, 1 / bheight * height)
    trafo.translate(x, y)
    tp = tp.transformed(trafo)
    return PathPatch(tp, facecolor=color, edgecolor=edgecolor)


def plot_motif(ax: mpl.axes.Axes,
               obj,  # list of weighted text
               graph_coords: Optional[Union[Dict, np.ndarray]] = None,
               width: float = 0.8,
               colors=None,
               width_per_character: float = 3.5,
               theme: str = "blank",
               **kwargs):
    Theme.set_theme(ax, theme)
    data = obj.data
    if colors is None:
        colors = ['#008000', '#cc0000', '#0000cc', '#ffb300']
        colors = {x: y for x, y in zip(["A", "T", "C", "G"], colors)}

    region = obj.region
    if graph_coords is None:
        graph_coords = init_graph_coords(region)

    # 在原始坐标轴上画motif
    if not data:
        logger.info("there is no any motif information to plot")
        return
    ymin, ymax, xmin, xmax = 0, 0, \
                             graph_coords[min(data.keys()) - region.start], \
                             graph_coords[max(data.keys()) - region.start] + (1 + width) / 2

    # 在放大区画motif
    axin_width = width_per_character * (xmax - xmin) / len(graph_coords)
    bbox_to_left = (xmax * 2 - xmin) / len(graph_coords)
    draw_left = False
    if axin_width + bbox_to_left > 1:
        bbox_to_left = max(xmin / len(graph_coords) - axin_width - .1, 0)
        draw_left = True

    # update insert_axes to matplotlib 3.7.2
    axins = ax.inset_axes(bounds=[bbox_to_left, 0, axin_width, 1], transform=ax.transAxes)

    start_site = min(list(data.keys()))
    site = 0
    for idx, vals in data.items():
        site = idx - start_site
        init_height_pos, init_height_neg = 0, 0
        for text, height in vals.items():
            text_shape = make_text_elements(text,
                                            x=site + (1 - width) / 2,
                                            y=init_height_neg if height < 0 else init_height_pos,
                                            width=width, height=height,
                                            color=colors.get(text, "blue"),
                                            edgecolor=colors.get(text, "blue"))

            if height > 0:
                init_height_pos += height
            else:
                init_height_neg += height
            ymin, ymax = min(ymin, init_height_neg), max(ymax, init_height_pos)
            axins.add_patch(text_shape)
    axins.set_xlim(0, site + (1 + width) / 2)
    axins.set_ylim(ymin, ymax)
    axins.tick_params(bottom=False, top=False, left=False, right=False)
    axins.set_xticklabels([])
    axins.set_yticklabels([])

    # draw box
    y_height = ymax - ymin
    y_center = (ymax + ymin) / 2
    y_scale = .25
    ax.plot(
        [xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymin + y_height * y_scale, ymin + y_height * y_scale, ymin],
        "black")

    ax.set_xlim(min(graph_coords), max(graph_coords))
    ax.set_ylim(ymin, ymax)

    # 画两条线
    if draw_left:
        pos = len(graph_coords) * ((bbox_to_left + axin_width) * 100 + 5) / 100
        x = xmin
    else:
        pos = len(graph_coords) * bbox_to_left
        x = xmax
    ax.arrow(x, ymin + y_height * y_scale / 2,
             pos - x, y_center - (ymin + y_height * y_scale),
             head_width=1, fc='k', ec='k', rasterized=True)


if __name__ == '__main__':
    pass

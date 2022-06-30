#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
This script contains the functions to draw different images
"""
import math
from copy import deepcopy
from typing import Dict, Optional, List, Union

import matplotlib as mpl
import numpy as np
import seaborn as sns
from matplotlib import pylab
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import gaussian_kde, zscore

from conf.config import DISTANCE_METRIC, CLUSTERING_METHOD
from conf.logger import logger
from sashimi.anno.theme import Theme
from sashimi.base.ReadSegments import ReadSegment
from sashimi.base.Stroke import Stroke
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.ReadDepth import ReadDepth
from sashimi.file.File import File
from sashimi.file.Reference import Reference


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
    return p0 * (1 - t) ** 3 + 3 * t * p1 * (1 - t) ** 2 + \
           3 * t ** 2 * (1 - t) * p2 + t ** 3 * p3


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
                      intron_scale=.5) -> np.array:
    u"""
    init the default
    :param region: the plot region
    :param exons: list of start and end sites of exons [[start, end], [start, end]]
    :param exon_scale: the scale of exon, default set to 1
    :param intron_scale: the scale of intron, default set to 0.5 -> the intron showed in final will be half size of real
    """
    graph_coords = np.zeros(len(region), dtype=int)

    if exons:
        # if there have exons, init graph_coords by exon and intron scales
        for i in range(0, exons[0][0] - region.start):
            graph_coords[i] = (i - 0) * intron_scale

        exons = __merge_exons__(exons)
        for i in range(0, len(exons)):
            exon = exons[i]

            if i > 0:
                intron = [exons[i - 1][1], exons[i][0]]

                for j in range(intron[0], intron[1]):
                    graph_coords[j - region.start] = graph_coords[intron[0] - region.start - 1] + (
                            j - intron[0] + 1) * intron_scale

            for j in range(exon[0], exon[1] + 1):
                graph_coords[j - region.start] = graph_coords[exon[0] - region.start - 1] + (
                        j - exon[0] + 1) * exon_scale

        intron = [exons[-1][-1], region.end]
        for i in range(intron[0], intron[1]):
            graph_coords[i - region.start] = graph_coords[intron[0] - region.start - 1] + (
                    i - intron[0] + 1) * intron_scale
    else:
        # if there is not any exons, just init graph_coords by region
        for i, j in enumerate(range(region.start, region.end + 1)):
            graph_coords[j - region.start] = i

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
    x_label = 'Genomic coordinate (%s), "%s" strand' % (region.chromosome, region.strand)

    if log_trans:
        if log_trans in ["2", "10"]:
            log_trans = f"log{log_trans}"
        x_label = f"{x_label}, y axis is {log_trans} transformed"

    ax.hlines(y=0, xmin=0, xmax=max(graph_coords), color="black", lw=1)
    ax.text(x=len(graph_coords) / 2, y=-2.8, s=x_label, fontsize=font_size, ha="center", va="top")

    bk = 1
    if not sequence:
        bk = len(graph_coords) // nx_ticks

    line_space = {}
    for i in range(0, len(graph_coords), bk):
        line_space[graph_coords[i]] = i

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
    ax.set_xlim(0, max(graph_coords))


def set_y_ticks(
        ax: mpl.axes.Axes,
        label: str,
        graph_coords: Union[Dict, np.array],
        max_used_y_val: Union[int, float],
        distance_between_label_axis: float = .1,
        n_y_ticks: int = 4,
        theme: str = "ticks",
        font_size: int = 5,
        show_y_label: bool = True,
        set_label_only: bool = False
):
    u"""
    The y ticks are formatted here
    @2019.03.31 add little check here to make sure the y axis shows the real value
    """
    # set y ticks, y label and label
    Theme.set_theme(ax, theme)

    if not set_label_only:
        ax.set_xlim(0, max(graph_coords))
        ax.set_ylim(- 0.5 * max_used_y_val, 1.2 * max_used_y_val)
        ax.spines["left"].set_bounds(0, max_used_y_val)

        universal_y_ticks = pylab.linspace(0, max_used_y_val, n_y_ticks + 1)

        curr_y_tick_labels = []
        for lab in universal_y_ticks:
            if lab <= 0:
                # Exclude label for 0
                curr_y_tick_labels.append("")
            else:
                curr_y_tick_labels.append("%.1f" % lab if lab % 1 != 0 else "%d" % lab)

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
        ax.text(
            x=-1 * distance_between_label_axis * max(graph_coords),
            y=max_used_y_val // 2, s=label,
            fontsize=font_size,
            ha="right"
        )


def set_focus(
        ax: mpl.axes.Axes,
        graph_coords: Union[Dict, np.array],
        focus: Dict[int, int]
):
    for left, right in focus.items():
        left, right = graph_coords[left], graph_coords[right]
        try:
            fill_x = [left, right, right, left]

            y1, y2 = ax.get_ylim()
            fill_y = [y1, y1, y2, y2]
            ax.fill(fill_x, fill_y, alpha=0.1, color='grey')
        except IndexError as err:
            logger.warning("focus region is out of bound: " + str(err))


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
            logger.warning("Indicator line is out of bound: " + str(err))


def plot_stroke(
        ax: mpl.axes.Axes,
        data: List[Stroke],
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        font_size: int = 5,
        distance_between_label_axis: Union[int, float] = .1,
        theme: str = "blank",
        **kwargs
):
    u"""
    plot stoke
    """

    strokes = sorted(data, key=lambda x: [x.start, x.end])

    for i, stroke in enumerate(strokes):
        ax.hlines(
            y=i,
            xmin=graph_coords[stroke.start],
            xmax=graph_coords[stroke.end],
            color=stroke.color, lw=2)
        ax.text(
            -1 * distance_between_label_axis * max(graph_coords),
            i + .2,
            stroke.label,
            fontsize=font_size,
            color=stroke.color
        )

    Theme.set_theme(ax, theme)
    ax.set_xlim(left=0, right=max(graph_coords))
    ax.set_ylim(bottom=-1, top=len(strokes))


def plot_reference(
        ax: mpl.axes.Axes,
        obj: Reference,
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        font_size: int = 5,
        show_gene: bool = False,
        show_id: bool = False,
        transcripts: Optional[List[str]] = None,
        remove_empty_transcripts: bool = False,
        color: Optional[str] = None,
        reverse_minus: bool = False,
        theme: str = "blank",
        y_loc: int = 0,
        exon_width: float = .3,
        plot_domain: bool = False,
        show_exon_id: bool = False,
        raster: bool = True
):
    u"""
    Plot the structure of reference
    :param ax:
    :param obj:
    :param graph_coords:
    :param font_size:
    :param show_gene:
    :param show_id:
    :param transcripts:
    :param remove_empty_transcripts:
    :param color:
    :param reverse_minus:
    :param theme:
    :param y_loc:
    :param exon_width:
    :param  plot_domain:
    :param show_exon_id:
    :param raster:
    :return:
    """
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

    for transcript in obj.data:
        # ignore the unwanted transcript
        if transcripts and not (set(transcripts) & transcript.ids()):
            continue

        # ignore transcripts without any exons
        if remove_empty_transcripts and not transcript.exons:
            continue

        strand = transcript.strand
        # @2018.12.20 add transcript id, based on fixed coordinates
        if transcript.transcript:
            if show_gene and transcript.gene and transcript.gene_id != transcript.transcript_id:
                if show_id:
                    ax.text( x=-1, y=y_loc + 0.25, s=transcript.gene_id, fontsize=font_size, ha="right")
                    ax.text(x=-1, y=y_loc - 0.25, s=transcript.transcript_id, fontsize=font_size, ha="right")
                else:
                    ax.text(x=-1, y=y_loc, s=transcript.gene + " | " + transcript.transcript,
                            fontsize=font_size, ha="right")
            else:
                ax.text(x=-1, y=y_loc - 0.1, s=transcript.transcript, fontsize=font_size, ha="right")

        # @2018.12.19
        # s and e is the start and end site of single exon
        # @2022.05.13
        # add index to avoid label overlapping of neighbor exon
        for ind, exon in enumerate(transcript.exons):
            s, e, strand = region.relative(exon.start), region.relative(exon.end), exon.strand
            x = [
                graph_coords[s], graph_coords[e],
                graph_coords[e], graph_coords[s]
            ]
            y = [
                y_loc - exon_width / 2, y_loc - exon_width / 2,
                y_loc + exon_width / 2, y_loc + exon_width / 2
            ]
            ax.fill(x, y, color, lw=.5, zorder=20)

            if show_exon_id:
                y_loc_offset = 0.1 if ind % 2 == 0 else - 0.2
                ax.text(x=(graph_coords[s] + graph_coords[s]) / 2, y=y_loc + y_loc_offset,
                        s=exon.name, fontsize=font_size / 2, ha="right")

        # @2018.12.21
        # change the intron range
        # Draw intron.
        # if transcript's category is interval, then don't plot the intron.
        if transcript.category != "interval":
            intron_sites = [
                graph_coords[region.relative(transcript.start)],
                graph_coords[region.relative(transcript.end)]
            ]
            ax.plot(intron_sites, [y_loc, y_loc], color=color, lw=0.5, rasterized=raster)

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
                    if strand == '+' or reverse_minus:
                        x = [loc - spread, loc, loc - spread]
                    else:
                        x = [loc + spread, loc, loc + spread]
                    y = [y_loc - exon_width / 5, y_loc, y_loc + exon_width / 5]
                    ax.plot(x, y, lw=.5, color=color, rasterized=raster)

        y_loc += 1  # if transcript.transcript else .5

        # here is plot domain
        # print(plot_domain)
        # print(obj.add_domain)
        # print(obj.domain)
        # print(transcript.transcript_id)

        if plot_domain and obj.domain and transcript.transcript_id in obj.domain.pep:
            current_domains = obj.domain.pep[transcript.transcript_id]

            for sub_current_domain in current_domains:
                if not (sub_current_domain.start <= region.end and
                        sub_current_domain.end >= region.start):
                    continue

                for sub_exon in sub_current_domain.exons:

                    for exon in sub_exon:
                        s, e = region.relative(exon.start), region.relative(exon.end)
                        if e < 0 or s > len(region):
                            continue
                        s = 0 if s < 0 else s
                        e = len(region) if e >= len(region) else e

                        x = [
                            graph_coords[s], graph_coords[e],
                            graph_coords[e], graph_coords[s]
                        ]
                        y = [
                            y_loc - exon_width / 4, y_loc - exon_width / 4,
                            y_loc + exon_width / 4, y_loc + exon_width / 4
                        ]
                        ax.fill(x, y, color, lw=.5, zorder=20, rasterized=raster)

                    # @2022.05.13
                    intron_relative_s = region.relative(min(map(lambda x: x.end, sub_exon)))
                    intron_relative_s = intron_relative_s if intron_relative_s >= 0 else 0

                    intron_relative_e = region.relative(max(map(lambda x: x.start, sub_exon)))
                    intron_relative_e = len(region) if intron_relative_e >= len(region) else intron_relative_e

                    intron_sites = [graph_coords[intron_relative_s], graph_coords[intron_relative_e]]
                    if len(sub_exon) != 1:
                        ax.plot(intron_sites, [y_loc, y_loc], color=color, lw=0.2, rasterized=raster)

                ax.text(x=-1, y=y_loc - 0.125, s=f"{sub_current_domain.gene}\n{transcript.transcript_id}",
                        fontsize=font_size / 2, ha="right")

                y_loc += 0.25
        # offset for next term.
        y_loc += 0.75

    # @2022.05.13
    # Set y lim using y_loc value.
    Theme.set_theme(ax, theme)
    ax.set_ylim(-.5, y_loc + .5)
    ax.set_xlim(0, max(graph_coords))


def plot_density(
        ax: mpl.axes.Axes,
        obj: Optional[File] = None,
        data: Optional[ReadDepth] = None,
        region: Optional[GenomicLoci] = None,
        graph_coords: Optional[Union[Dict, np.ndarray]] = None,
        color="blue",
        font_size: int = 8,
        show_junction_number: bool = True,
        junction_number_font_size: int = 5,
        n_y_ticks: int = 4,
        distance_between_label_axis: float = .1,
        show_y_label: bool = True,
        y_label: str = "",
        theme: str = "ticks_blank",
        max_used_y_val: Optional[float] = None,
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
    :param max_used_y_val: used to set same max y axis
    :return:
    """
    if obj:
        assert obj.region is not None, "please load data first"
        region = obj.region
        if graph_coords is None:
            graph_coords = init_graph_coords(region)

        data = deepcopy(obj.data)
        data.transform(obj.log_trans)

        if not y_label:
            y_label = obj.label

    elif not (region and data):
        raise ValueError("please input obj or region and data")
    wiggle = data.wiggle

    if max_used_y_val is None:
        max_used_y_val = max(wiggle)
        if max_used_y_val % 2 == 1:
            max_used_y_val += 1

    jxns = data.junctions_dict

    y_max = max_used_y_val
    y_min = -.5 * y_max

    # Reduce memory footprint by using incremented graphcoords.
    compressed_x = []
    compressed_wiggle = []

    u"""
    @2019.01.04
    If there is no bam file, use half of y axis as the upper bound of exon 
    And draw a white point to maintain the height of the y axis
    """
    for i in range(len(graph_coords)):
        compressed_wiggle.append(wiggle[i])
        compressed_x.append(graph_coords[i])

    ax.fill_between(compressed_x, compressed_wiggle, y2=0, color=color, lw=0, step="post")

    if jxns:
        # sort the junctions by intron length for better plotting look
        jxns_sorted_list = sorted(jxns.keys(), key=lambda x: x.end - x.start, reverse=True)

        if not jxns:
            max_junction_count, min_junction_count = 0, 0
        else:
            max_junction_count = max(jxns.values())
            min_junction_count = min(jxns.values())
        junction_count_gap = max_junction_count - min_junction_count

        current_height = -3 * y_min / 4

        for plotted_count, jxn in enumerate(jxns_sorted_list):
            leftss, rightss = jxn.start, jxn.end

            # @2018.12.19
            # set junctions coordinate here
            # the junction out of boundaries, set the boundaries as coordinate
            ss1_idx, ss1_modified = get_limited_index(leftss - region.start, len(graph_coords))
            ss2_idx, ss2_modified = get_limited_index(rightss - region.start, len(graph_coords))

            u"""
            @2019.01.14
            add two new variables to make it clear which one is index, which one is genomic site 
            """
            ss1 = graph_coords[ss1_idx]
            ss2 = graph_coords[ss2_idx]

            # draw junction on bottom
            if plotted_count % 2 == 0:

                pts = [
                    (ss1, 0 if not ss1_modified else -current_height),
                    (ss1, -current_height),
                    (ss2, -current_height),
                    (ss2, 0 if not ss2_modified else -current_height)
                ]
                midpt = cubic_bezier(pts, .5)

            # draw junction on top
            else:

                left_dens = wiggle[ss1_idx]
                right_dens = wiggle[ss2_idx]

                """
                @2019.01.04
    
                If there is no bam, lower half of y axis as the height of junctions
                """
                pts = [
                    (ss1, left_dens if not ss1_modified else left_dens + current_height),
                    (ss1, left_dens + current_height),
                    (ss2, right_dens + current_height),
                    (ss2, right_dens if not ss2_modified else right_dens + current_height)
                ]

                midpt = cubic_bezier(pts, .5)

            if show_junction_number:
                t = ax.text(
                    midpt[0], midpt[1],
                    '{0}'.format(round(jxns[jxn], 2)),
                    fontsize=junction_number_font_size,
                    ha='center', va='center',
                    backgroundcolor='w'
                )

                # @2018.12.19 transparent background
                t.set_bbox(dict(alpha=0))

            a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])

            """
            @2018.12.26
            scale the junctions line width
            """
            if junction_count_gap > 0:
                line_width = (jxns[jxn] - min_junction_count) / junction_count_gap
            else:
                line_width = 0

            ax.add_patch(PathPatch(a, ec=color, lw=line_width + 0.2, fc='none'))

    if obj and obj.title:
        ax.text(
            max(graph_coords) - len(obj.title),
            max_used_y_val,
            obj.title,
            color=color,
            fontsize=font_size
        )

    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=max_used_y_val,
        n_y_ticks=n_y_ticks,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
    )


def plot_side_plot(
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

    data = deepcopy(obj.data)
    data.transform(obj.log_trans)

    if not y_label:
        y_label = obj.label

    if not show_y_label:
        y_label = ""

    plus, minus = data.plus, data.minus

    max_height = max(plus)
    min_height = min(minus)
    max_val = max(max_height, abs(min_height))

    for label, array_plot in zip(['plus', 'minus'], [plus, minus]):
        if strand_choice != "all" and label != strand_choice:
            continue

        array_hist = np.repeat(graph_coords, np.abs(array_plot).astype(np.int))
        try:
            kde = gaussian_kde(array_hist)
            fit_value = kde.pdf(graph_coords)
        except ValueError:
            # logger.warning(err)
            # logger.warning(traceback.format_exc())
            continue

        fit_value = fit_value / fit_value.max()
        if label == 'plus':
            ax.plot(graph_coords, fit_value * array_plot.max(), c=color, lw=1)
            ax.bar(range(len(graph_coords)), array_plot, color=color, rasterized=raster)
        else:
            ax.plot(graph_coords, fit_value * array_plot.min(), c=color, lw=1)
            ax.bar(range(len(graph_coords)), array_plot, color=color, rasterized=raster)

    # set the y limit
    # set y ticks, y label and label
    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=1.1 * max_val,
        n_y_ticks=n_y_ticks,
        font_size=font_size,
        show_y_label=False,
        distance_between_label_axis=distance_between_label_axis
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
        raster: bool = True
):
    u"""

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
    """
    labels = list(data.keys())
    mtx = np.array([x.wiggle for x in data.values()])

    if clustering and len(mtx) > 1:
        assert clustering_method in CLUSTERING_METHOD, f"clustering_method {clustering_method} is not supported."
        assert distance_metric in DISTANCE_METRIC, f"distance_metric {distance_metric} is not supported"

        order = dendrogram(
            linkage(mtx, method=clustering_method, metric=distance_metric),
            orientation='right'
        )

        mtx = mtx[order["leaves"], :]
        labels = [labels[x] for x in order["leaves"]]

    if do_scale:
        """
        y = (x – mean) / standard_deviation
        """
        mtx = zscore(mtx, axis=1)

    if not show_y_label:
        labels = False

    sns.heatmap(
        mtx,
        ax=ax, cmap=color,
        cbar_ax=cbar_ax,
        xticklabels=False,
        yticklabels=labels,
        center=True,
        rasterized=raster
    )

    ax.tick_params(axis='both', which='major', labelsize=font_size)
    cbar_ax.tick_params(labelsize=font_size)

    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=mtx.shape[0],
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        set_label_only=True
    )
    if raster:
        ax.set_rasterization_zorder(0)


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
    max_x_val = 0

    for ylab, val in data.items():
        attr = line_attrs.get(ylab, {}) if line_attrs else {}

        x, y = [], []
        for i, j in enumerate(val.wiggle):
            x.append(graph_coords[i])
            y.append(j)
        ax.plot(x, y, label=ylab if show_y_label else "", **attr)

        max_y_val = max(max_y_val, max(val.wiggle))
        max_x_val = max(max_x_val, len(val.wiggle))

    ax.tick_params(axis='both', which='major', labelsize=font_size)
    ax.tick_params(labelsize=font_size)

    if show_legend:
        ax.legend(loc=legend_position,
                  ncol=int(len(data) / 1.5) if legend_ncol <= 0 else legend_ncol,
                  fancybox=False, shadow=False)

    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=max_y_val if max_used_y_val is None else max_used_y_val,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        n_y_ticks=n_y_ticks
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
        raster: bool = False
):
    u"""

    :param show_y_label:
    :param n_y_ticks:
    :param feature_color:
    :param distance_between_label_axis:
    :param ax:
    :param obj:
    :param graph_coords:
    :param y_label:
    :param exon_color:
    :param intron_color:
    :param exon_width:
    :param font_size:
    :param theme:
    :param raster:
    :return:
    """

    assert len(obj) == 1, "IGV-like plot only support one file"

    obj = list(obj.values())[0]
    assert obj.region is not None, "please load data first"

    region = obj.region
    if graph_coords is None:
        graph_coords = init_graph_coords(region)

    if not y_label:
        y_label = obj.label

    y_loc = 1.0

    for c_ind_list in obj.get_index():
        for c_ind in c_ind_list:
            c_data = obj.data[c_ind]
            # skip truncated reads in the given region
            if c_data.start < region.start or c_data.end > region.end:
                continue

            for exon in c_data.exons:
                s, e, strand = region.relative(exon.start), region.relative(exon.end), exon.strand
                if e < 0 or s > len(region):
                    continue

                s = 0 if s < 0 else s
                e = len(region) - 1 if e >= len(region) else e

                x = [
                    graph_coords[s], graph_coords[e],
                    graph_coords[e], graph_coords[s]
                ]
                y = [
                    y_loc - exon_width, y_loc - exon_width,
                    y_loc + exon_width, y_loc + exon_width
                ]

                ax.fill(x, y, 'k' if not exon_color else exon_color, lw=.5, zorder=20)

            for intron in c_data.introns:
                s, e, strand = region.relative(intron.start), region.relative(intron.end), intron.strand

                if e < 0 or s > len(region):
                    continue
                s = 0 if s < 0 else s
                e = len(region) - 1 if e >= len(region) else e
                intron_sites = [
                    graph_coords[s],
                    graph_coords[e]
                ]
                ax.plot(intron_sites, [y_loc, y_loc],
                        color="#4d4d4d" if not intron_color else intron_color,
                        lw=0.2, rasterized=raster)

            for feature in c_data.features:
                if feature.start == -1:
                    continue

                s, e, strand = region.relative(feature.start), region.relative(feature.end), feature.strand
                is_site = False

                if s == e:
                    s = s - 1
                    is_site = True

                if e < 0 or s > len(region):
                    continue

                s = 0 if s < 0 else s
                e = len(region) - 1 if e >= len(region) else e

                x = [
                    graph_coords[s], graph_coords[e],
                    graph_coords[e], graph_coords[s]
                ]
                width_ratio = 0.75 if not is_site else 1.5

                y = [
                    y_loc - exon_width * width_ratio, y_loc - exon_width * width_ratio,
                    y_loc + exon_width * width_ratio, y_loc + exon_width * width_ratio
                ]
                ax.fill(x, y, 'r' if not feature_color else feature_color, lw=.2, zorder=20)
                if is_site:
                    ax.scatter(graph_coords[s],
                               y_loc + exon_width * width_ratio,
                               c='b' if not feature_color else feature_color,
                               s=.5,
                               linewidths=(0,),
                               rasterized=raster
                               )

        y_loc += 1
        #

    set_y_ticks(
        ax, label=y_label, theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=y_loc,
        n_y_ticks=n_y_ticks,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
    )


if __name__ == '__main__':
    from matplotlib import pyplot as plt


    def test_density():
        from sashimi.file.Bam import Bam
        fig, ax = plt.subplots()
        region = GenomicLoci("chr1", 1270656, 1284730, "+")
        bam = Bam.create("../example/bams/1.bam")
        bam.load(region, log_trans="2")
        plot_density(ax, bam)
        plt.savefig("plot_density_bam.png")


    def test_bw():
        from sashimi.file.Bigwig import Bigwig

        bw = Bigwig.create("../example/bws/1.bw", title="test")
        bw.load(GenomicLoci("chr1", 1270656, 1284730, "+"))
        fig, ax = plt.subplots()
        plot_density(ax, bw)
        plt.savefig("plot_density_bw.png")


    def test_ref():
        region = GenomicLoci("chr1", 1270656, 1284730, "+")
        fig, ax = plt.subplots()
        ref = Reference.create("../example/example.sorted.gtf.gz")
        ref.load(region, domain=True)

        ref.add_interval(
            interval_file="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.bed.gz",
            interval_label="PolyASite"
        )

        ref.add_interval(
            interval_file="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.simple.bed.gz",
            interval_label="PolyASite_simple"
        )

        plot_reference(ax, ref,
                       show_gene=True,
                       show_id=True,
                       plot_domain=True,
                       show_exon_id=True
                       )

        plt.savefig("plot_reference.pdf")


    def test_depth():
        from sashimi.file.Depth import Depth
        region = GenomicLoci("chr1", 1270656, 1284730, "+")
        depth = Depth.create("../example/depth.bgz")
        depth.load(region)
        fig, ax = plt.subplots(nrows=len(depth))
        idx = 0
        for x, y in depth.items():
            plot_density(ax[idx], region=region, data=y, y_label=x,
                         theme="ticks_blank" if idx < len(depth) - 1 else "ticks")
            idx += 1

        plt.savefig("plot_density_depth.png")


    def test_heatmap_and_line():
        from matplotlib import gridspec
        from sashimi.file.Bam import Bam

        region = GenomicLoci("chr1", 1270656, 1284730, "+")
        graph_coords = init_graph_coords(region)
        data = {}
        attrs = {}
        colors = ["red", "blue", "yellow", "green", "grey"]
        for i in range(1, 5):
            bam = Bam.create(f"../example/bams/{i}.bam")
            bam.load(region)
            data[str(i)] = bam.data
            attrs[str(i)] = {"color": colors[i]}

        gs = gridspec.GridSpec(1, 2, width_ratios=(.99, .01), wspace=0.01, hspace=.15)
        ax_var = plt.subplot(gs[:, 0]),
        cbar_ax = plt.subplot(gs[:, 1])
        plot_heatmap(ax_var[0], cbar_ax, data, do_scale=True, clustering=True)
        plt.savefig("plot_heatmap.png")

        fig, ax = plt.subplots()
        plot_line(ax, data, line_attrs=attrs, graph_coords=graph_coords)
        plt.savefig("plot_line.png")


    def test_graph_coord():
        from matplotlib import gridspec
        region = GenomicLoci(chromosome="1", start=100, end=800, strand="+")
        exon_width = .5
        exons = [
            [150, 300], [320, 450],
            [460, 500], [470, 500],
            [470, 600], [700, 800],
        ]

        gs = gridspec.GridSpec(2, 1)
        ax = plt.subplot(gs[0, 0])
        for ind, exon in enumerate(exons):
            s, e, strand = exon[0], exon[1], "+"
            x = [s, e, e, s]
            x = [i - region.start for i in x]
            y = [
                ind - exon_width / 2, ind - exon_width / 2,
                ind + exon_width / 2, ind + exon_width / 2
            ]
            ax.fill(x, y, 'k', lw=.5, zorder=20)

        ax.set_xbound(0, len(region))
        graph_coords = init_graph_coords(region, exons, intron_scale=.5)

        ax = plt.subplot(gs[1, 0])
        for ind, exon in enumerate(exons):
            s, e, strand = exon[0] - region.start, exon[1] - region.start, "+"
            x = [graph_coords[s], graph_coords[e], graph_coords[e], graph_coords[s]]
            y = [
                ind - exon_width / 2, ind - exon_width / 2,
                ind + exon_width / 2, ind + exon_width / 2
            ]
            ax.fill(x, y, 'k', lw=.5, zorder=20)
        ax.set_xbound(0, len(region))
        plt.savefig("plot_graph_coords.png")


    def test_set_x_ticks():
        region = GenomicLoci(chromosome="1", start=100, end=120, strand="+")
        fig, ax = plt.subplots()

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(left=False)
        set_x_ticks(ax, region, font_size=10, sequence={110: "A", 106: "C"})
        plt.savefig("plot_x_ticks.png")


    def test_igv_plot():
        from sashimi.base.ReadSegments import ReadSegment
        fig, ax = plt.subplots()
        region = GenomicLoci("chr1", 13362, 29900, "+")
        rs = ReadSegment.create("../example/bams/WASH7P.bam")
        rs.load(region, features={"m6a": "ma", "real_strand": "rs", "polya": "pa"})
        plot_igv_like(ax, rs, y_label="fl")
        plt.savefig("test_igv_plot.pdf")


    def test_igv_plot2():
        from sashimi.base.ReadSegments import ReadSegment
        fig, ax = plt.subplots()
        region = GenomicLoci("chr1", 1270656, 1284730, "+")
        rs = ReadSegment.create("../example/bams/0.bam")
        rs.load(region, features={"m6a": "ma", "real_strand": "rs", "polya": "pa"})
        plot_igv_like(ax, rs, y_label="fl")
        plt.savefig("test_igv_plot.2.pdf")


    test_igv_plot()
    test_igv_plot2()

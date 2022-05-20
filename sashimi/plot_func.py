#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
This script contains the functions to draw different images
"""
import math
from copy import deepcopy
from typing import Dict, Optional, List

import matplotlib as mpl
import numpy as np
import seaborn as sns

from matplotlib import pylab
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import zscore

from sashimi.anno.theme import Theme
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.ReadDepth import ReadDepth
from sashimi.file.File import File
from sashimi.file.Reference import Reference

from conf.heatmap import DISTANCE_METRIC, CLUSTERING_METHOD


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


# def set_x_ticks(read_depth_object, ax_var, graph_coords, chromosome, strand, logtrans = None, nx_ticks = 4, font_size = 4):
#     ax_var.xaxis.set_ticks_position('bottom')
#
#     # @2018.12.19 unnecessary text in figure
#
#     xlabel = 'Genomic coordinate (%s), "%s" strand' % (chromosome, strand)
#
#     if logtrans in (2, 10):
#         xlabel = xlabel + ", y axis is log%d transformed" % logtrans
#
#     pylab.xlabel(xlabel, fontsize=font_size)
#
#     bk = 1
#     if not read_depth_object.sequence:
#         bk = len(graph_coords) // nx_ticks
#
#     linspace, ticks = [], []
#     for i in range(0, len(graph_coords), bk):
#         linspace.append(graph_coords[i])
#
#         temp_txs = read_depth_object.start + i
#         if read_depth_object.sequence:
#             if (i - read_depth_object.start) % nx_ticks == 0:
#                 temp_txs = "{}\n{}".format(read_depth_object.start + i, read_depth_object.sequence[i])
#             else:
#                 temp_txs = "\n{}".format(read_depth_object.sequence[i])
#         ticks.append(temp_txs)
#
#     # ax_var.spines['right'].set_color('none')
#     # ax_var.spines['top'].set_color('none')
#     # pylab.yticks([])
#     ax_var.spines['bottom'].set_color('black')
#     ax_var.xaxis.set_ticks_position('bottom')
#     pylab.xticks(linspace, ticks, fontsize=font_size)


def set_y_ticks(
        ax_var: mpl.axes.Axes,
        label: str,
        universal_y_ticks,
        distance_between_label_axis: float,
        font_size=5,
        show_y_label: bool = True
):
    u"""
    The y ticks are formatted here
    @2019.03.31 add little check here to make sure the y axis shows the real value
    """
    """
    Plot y labels
    @2018.12.20 using BAM label as ylabel
    @2019.01.04 change the standards of distance between ylabel and y-axis
    """
    if show_y_label:
        ax_var.set_ylabel(
            label,
            fontsize=font_size,
            va="center",
            labelpad=distance_between_label_axis,  # the distance between ylabel with axis
            rotation="horizontal"
        )

    curr_y_tick_labels = []
    for label in universal_y_ticks:
        if label <= 0:
            # Exclude label for 0
            curr_y_tick_labels.append("")
        else:
            curr_y_tick_labels.append("%.1f" % label if label % 1 != 0 else "%d" % label)

    u"""
    @2019.01.04
    If there is no bam file, draw a blank y-axis 
    """
    ax_var.set_yticks(universal_y_ticks)
    ax_var.set_yticklabels(curr_y_tick_labels, fontsize=font_size)
    ax_var.yaxis.set_ticks_position('left')


def init_graph_coords(region: GenomicLoci) -> np.array:
    graph_coords = np.zeros(len(region), dtype=int)
    for i, j in enumerate(range(region.start, region.end + 1)):
        graph_coords[j - region.start] = i
    return graph_coords


def plot_reference(
        ax: mpl.axes.Axes,
        obj: File,
        graph_coords: Optional[dict] = None,
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
        show_exon_id: bool = False
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
    :return:
    """
    Theme.set_theme(ax, theme)
    region = obj.region

    data = sorted(obj.data)
    if not data:
        return

    if not graph_coords:
        graph_coords = init_graph_coords(region)

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
                    ax.text(
                        x=-1, y=y_loc + 0.25, s=transcript.gene_id,
                        fontsize=font_size,
                        ha="right"
                    )

                    ax.text(
                        x=-1, y=y_loc - 0.25,
                        s=transcript.transcript_id,
                        fontsize=font_size,
                        ha="right"
                    )
                else:
                    ax.text(
                        x=-1, y=y_loc,
                        s=transcript.gene + " | " + transcript.transcript,
                        fontsize=font_size,
                        ha="right"
                    )
            else:
                ax.text(x=-1, y=y_loc - 0.1,
                        s=transcript.transcript,
                        fontsize=font_size,
                        ha="right")

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
            ax.fill(x, y, 'k' if not color else color, lw=.5, zorder=20)

            if show_exon_id:
                y_loc_offset = 0.1 if ind % 2 == 0 else - 0.2
                ax.text(x=(graph_coords[s] + graph_coords[s]) / 2,
                        y=y_loc + y_loc_offset,
                        s=exon.name,
                        fontsize=font_size / 2,
                        ha="right")

        # @2018.12.21
        # change the intron range
        # Draw intron.
        # if transcript's category is interval, then don't plot the intron.
        if transcript.category != "interval":
            intron_sites = [
                graph_coords[region.relative(transcript.start)],
                graph_coords[region.relative(transcript.end)]
            ]
            ax.plot(intron_sites, [y_loc, y_loc], color='k' if not color else color, lw=0.5)

            # @2018.12.23 fix intron arrows issues
            # Draw intron arrows.
            max_ = graph_coords[region.relative(transcript.end)]
            min_ = graph_coords[region.relative(transcript.start)]
            length = max_ - min_
            narrows = math.ceil(length / max(graph_coords) * 50)

            if narrows > 0:
                spread = .2 * length / narrows
                # print(f"{spread} = .2 * {length} / {narrows}")
                for i in range(narrows):
                    loc = float(i) * length / narrows + graph_coords[region.relative(transcript.start)]
                    if strand == '+' or reverse_minus:
                        x = [loc - spread, loc, loc - spread]
                    else:
                        x = [loc + spread, loc, loc + spread]
                    y = [y_loc - exon_width / 5, y_loc, y_loc + exon_width / 5]
                    ax.plot(x, y, lw=.5, color='k')

        y_loc += 1  # if transcript.transcript else .5

        # here is plot domain
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
                        ax.fill(x, y, 'k' if not color else color, lw=.5, zorder=20)

                    # print(sub_exon)
                    # @2022.05.13
                    #
                    intron_relative_s = region.relative(
                        min(map(lambda x: x.end, sub_exon))
                    )
                    intron_relative_s = intron_relative_s if intron_relative_s >= 0 else 0

                    intron_relative_e = region.relative(
                        max(map(lambda x: x.start, sub_exon))
                    )
                    intron_relative_e = len(region) if intron_relative_e >= len(region) else intron_relative_e

                    intron_sites = [
                        graph_coords[intron_relative_s],
                        graph_coords[intron_relative_e]
                    ]
                    if len(sub_exon) != 1:
                        ax.plot(intron_sites, [y_loc, y_loc], color='k' if not color else color, lw=0.2)

                ax.text(
                    x=-1, y=y_loc - 0.125, s=f"{sub_current_domain.gene}\n{transcript.transcript_id}",
                    fontsize=font_size/2,
                    ha="right"
                )

                y_loc += 0.25
        # offset for next term.
        y_loc += 0.75

    # @2022.05.13
    # Set ylim using y_loc value.
    ax.set_ylim(-.5, y_loc + .5)


def plot_density(
        ax: mpl.axes.Axes,
        obj: Optional[File] = None,
        data: Optional[ReadDepth] = None,
        region: Optional[GenomicLoci] = None,
        graph_coords: Optional[dict] = None,
        color="blue",
        font_size: int = 8,
        show_junction_number: bool = True,
        junction_number_font_size: int = 5,
        ny_ticks: int = 4,
        distance_between_label_axis: float = .3,
        show_y_label: bool = True,
        y_label: str = "",
        theme: str = "ticks_blank"
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
    :param ny_ticks: number of y ticks
    :param junction_number_font_size:
    :param obj: Bam or Bigwig object
    :param color: color for this density plot
    :param show_y_label: whether to show y-axis label
    :param y_label: the text of y-axis title
    :param theme: the theme name
    :return:
    """
    if obj:
        assert obj.region is not None, "please load data first"
        region = obj.region
        if not graph_coords:
            graph_coords = init_graph_coords(region)

        data = deepcopy(obj.data)
        data.transform(obj.log_trans)

        if not y_label:
            y_label = obj.label
    elif not (region and data):
        raise ValueError("please input obj or region and data")
    wiggle = data.wiggle

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

    # set y ticks, y label and label
    Theme.set_theme(ax, theme)
    ax.set_xbound(0, max(graph_coords))
    ax.set_ybound(lower=- 0.5 * max_used_y_val, upper=1.2 * max_used_y_val)
    ax.spines["left"].set_bounds(0, max_used_y_val)

    universal_y_ticks = pylab.linspace(0, max_used_y_val, ny_ticks + 1)
    set_y_ticks(
        ax,
        label=y_label,
        universal_y_ticks=universal_y_ticks,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
    )


def plot_heatmap(
        ax: mpl.axes.Axes,
        cbar_ax: mpl.axes.Axes,
        data: Dict[str, ReadDepth],
        color="viridis",
        font_size: int = 8,
        distance_between_label_axis: float = .3,
        show_y_label: bool = True,
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

    if y_label:
        ax.set_ylabel(
            y_label,
            fontsize=font_size,
            va="center",
            labelpad=distance_between_label_axis,  # the distance between y label with axis
            rotation="horizontal"
        )

    if raster:
        ax.set_rasterization_zorder(0)


def plot_line(
        ax: mpl.axes.Axes,
        data: Dict[str, ReadDepth],
        font_size: int = 8,
        distance_between_label_axis: float = .3,
        show_y_label: bool = True,
        y_label: str = "",
        line_attrs: Optional[Dict[str, Dict]] = None,
        theme: str = "ticks_blank",
        ny_ticks: int = 4,
        show_legend: bool = False
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
    :param line_attrs: dict contains attributes like color and type for each line
    :param theme: the theme name
    :param ny_ticks: the number of y ticks
    :param show_legend: whether to show legend
    """
    max_used_y_val = 0
    max_used_x_val = 0
    for ylab, val in data.items():
        attr = line_attrs.get(ylab, {}) if line_attrs else {}
        ax.plot(range(len(val.wiggle)), val.wiggle, label=ylab if show_y_label else "", **attr)
        max_used_y_val = max(max_used_y_val, max(val.wiggle))
        max_used_x_val = max(max_used_x_val, len(val.wiggle))

    ax.tick_params(axis='both', which='major', labelsize=font_size)
    ax.tick_params(labelsize=font_size)

    if show_legend:
        ax.legend()

    if y_label:
        ax.set_ylabel(
            y_label,
            fontsize=font_size,
            va="center",
            labelpad=distance_between_label_axis,  # the distance between y label with axis
            rotation="horizontal"
        )

    Theme.set_theme(ax, theme)

    ax.set_xbound(0, max_used_x_val)
    ax.set_ybound(lower=- 0.5 * max_used_y_val, upper=1.2 * max_used_y_val)
    ax.spines["left"].set_bounds(0, max_used_y_val)
    universal_y_ticks = pylab.linspace(0, max_used_y_val, ny_ticks + 1)
    set_y_ticks(
        ax,
        label=y_label,
        universal_y_ticks=universal_y_ticks,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
    )


if __name__ == '__main__':
    from matplotlib import pyplot as plt

    from sashimi.file.Bam import Bam
    from sashimi.file.Bigwig import Bigwig
    from sashimi.file.Depth import Depth

    # fig, ax = plt.subplots()
    #
    region = GenomicLoci("chr1", 1270656, 1284730, "+")
    # bam = Bam.create("../example/bams/1.bam")
    # bam.load(region, log_trans="2")
    # plot_density(ax, bam)
    # plt.savefig("plot_density_bam.png")
    #
    # bw = Bigwig.create("../example/bws/1.bw", title="test")
    # bw.load(GenomicLoci("chr1", 1270656, 1284730, "+"))
    #
    # plot_density(ax, bw)
    # plt.savefig("plot_density_bw.png")
    #
    # fig, ax = plt.subplots()
    # ref = Reference.create("../example/example.sorted.gtf.gz")
    # ref.load(region, domain=True)
    #
    # ref.add_interval(
    #     interval_file="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.bed.gz",
    #     interval_label="PolyASite"
    # )
    #
    # ref.add_interval(
    #     interval_file="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.simple.bed.gz",
    #     interval_label="PolyASite_simple"
    # )
    #
    # plot_reference(ax, ref,
    #                show_gene=True,
    #                show_id=True,
    #                plot_domain=True,
    #                show_exon_id=True
    #                )
    #
    # plt.savefig("plot_reference.pdf")

    # depth = Depth.create("../example/depth.bgz")
    # depth.load(region)
    # fig, ax = plt.subplots(nrows=len(depth))
    # idx = 0
    # for x, y in depth.items():
    #     plot_density(ax[idx], region=region, data=y, y_label=x, theme="ticks_blank" if idx < len(depth) - 1 else "ticks")
    #     idx += 1
    #
    # plt.savefig("plot_density_depth.png")

    data = {}
    attrs = {}
    colors = ["red", "blue", "yellow", "green", "grey"]
    for i in range(1, 5):
        bam = Bam.create(f"../example/bams/{i}.bam")
        bam.load(region)
        data[str(i)] = bam.data
        attrs[str(i)] = {"color": colors[i]}

    from matplotlib import gridspec
    # gs = gridspec.GridSpec(1, 2, width_ratios=(.99, .01), wspace=0.01, hspace=.15)
    # ax_var = plt.subplot(gs[:, 0]),
    # cbar_ax = plt.subplot(gs[:, 1])
    #
    # plot_heatmap(ax_var[0], cbar_ax, data, do_scale=True, clustering=True)
    # plt.savefig("plot_heatmap.png")

    fig, ax = plt.subplots()
    plot_line(ax, data, line_attrs=attrs)
    plt.savefig("plot_line.png")

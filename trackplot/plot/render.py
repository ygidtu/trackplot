#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
All plot rendering functions.

Dependencies: coord (axis helpers), limits (y-limit computation), utils
"""

import math
from typing import Dict, List, Optional, Union

import matplotlib as mpl
import numpy as np
import seaborn as sns
from adjustText import adjust_text
from loguru import logger
from matplotlib import pylab as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch, Polygon
from matplotlib.path import Path
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import gaussian_kde, zscore

from trackplot.anno.theme import Theme
from trackplot.base.GenomicLoci import GenomicLoci
from trackplot.base.ReadDepth import ReadDepth
from trackplot.base.Stroke import Stroke
from trackplot.conf.config import CLUSTERING_METHOD, DISTANCE_METRIC
from trackplot.file.Annotation import Annotation
from trackplot.file.File import File
from trackplot.file.HiCMatrixTrack import HiCTrack
from trackplot.file.ReadSegments import ReadSegment
from trackplot.plot.coord import init_graph_coords, set_y_ticks
from trackplot.plot.limits import _compute_y_limits
from trackplot.plot.utils import (cubic_bezier, get_limited_index,
                                  make_text_elements)

# ============================================================================
# Stroke
# ============================================================================


def plot_stroke(
    ax,
    data: List[Stroke],
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    font_size: int = 5,
    theme: str = "blank",
    **kwargs,
):
    strokes = sorted(data, key=lambda x: [x.start, x.end])

    for i, stroke in enumerate(strokes):
        try:
            ax.hlines(
                y=i,
                xmin=graph_coords[stroke.start],
                xmax=graph_coords[stroke.end],
                color=stroke.color,
                lw=2,
            )
            ax.text(
                -0.1 * max(graph_coords),
                i + 0.2,
                stroke.label,
                fontsize=font_size,
                color=stroke.color,
            )
        except IndexError:
            logger.info(
                f"stroke {stroke} is out of bound, this stroke won't show up in final plot"
            )

    Theme.set_theme(ax, theme)
    ax.set_xlim(min(graph_coords), max(graph_coords))
    ax.set_ylim(bottom=-1, top=len(strokes))


# ============================================================================
# Annotation
# ============================================================================


def plot_annotation(
    ax,
    obj: Annotation,
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    font_size: int = 5,
    show_gene: bool = False,
    show_id: bool = False,
    color: Optional[str] = None,
    theme: str = "blank",
    y_loc: int = 0,
    exon_width: float = 0.3,
    plot_domain: bool = False,
    show_exon_id: bool = False,
    raster: bool = True,
    **kwargs,
):
    Theme.set_theme(ax, theme)
    ax.set_xlim(min(graph_coords), max(graph_coords))

    region = obj.region
    data = sorted(obj.data)
    if not data:
        return

    if graph_coords is None:
        graph_coords = init_graph_coords(region)

    color = "k" if not color else color

    patches, arrows = [], {}
    for transcript in data:
        strand = transcript.strand
        if transcript.transcript or transcript.transcript_id:
            show_text = transcript.transcript
            if not transcript.transcript:
                logger.warning(
                    "there is not transcript_name, using transcript_id instead"
                )
                show_text = transcript.transcript_id

            if (
                show_gene
                and transcript.gene
                and transcript.gene_id != transcript.transcript_id
            ):
                if show_id:
                    ax.text(
                        x=-0.01 * max(graph_coords),
                        y=y_loc + 0.25,
                        s=transcript.gene_id,
                        fontsize=font_size,
                        ha="right",
                    )
                    ax.text(
                        x=-0.01 * max(graph_coords),
                        y=y_loc - 0.25,
                        s=transcript.transcript_id,
                        fontsize=font_size,
                        ha="right",
                    )
                else:
                    ax.text(
                        x=-0.01 * max(graph_coords),
                        y=y_loc,
                        s=transcript.gene + " | " + show_text,
                        fontsize=font_size,
                        ha="right",
                    )
            else:
                ax.text(
                    x=-1, y=y_loc - 0.1, s=show_text, fontsize=font_size, ha="right"
                )
        else:
            logger.warning("there is no transcript_name and transcript_id")

        for ind, exon in enumerate(transcript.exons):
            s, e = region.relative(exon.start), region.relative(exon.end)
            if s < 0 or e < 0:
                continue

            patches.append(
                plt.Rectangle(
                    (graph_coords[s], y_loc - exon_width / 2),
                    width=graph_coords[e] - graph_coords[s],
                    height=exon_width,
                    lw=0.5,
                    zorder=20,
                    rasterized=raster,
                    color=color,
                )
            )
            if show_exon_id and exon.name:
                y_loc_offset = 0.1 if ind % 2 == 0 else -0.2
                ax.text(
                    x=(graph_coords[s] + graph_coords[e]) / 2,
                    y=y_loc + y_loc_offset,
                    s=exon.name,
                    fontsize=font_size / 2,
                    ha="center",
                    color="red",
                )

        if transcript.category != "interval":
            intron_sites = [
                graph_coords[region.relative(transcript.start)],
                graph_coords[region.relative(transcript.end)],
            ]
            ax.plot(intron_sites, [y_loc, y_loc], color=color, lw=0.5)

            max_ = graph_coords[region.relative(transcript.end)]
            min_ = graph_coords[region.relative(transcript.start)]
            length = max_ - min_
            narrows = math.ceil(length / max(graph_coords) * 50)
            if narrows > 0:
                spread = 0.2 * length / narrows
                for i in range(narrows):
                    loc = (
                        float(i) * length / narrows
                        + graph_coords[region.relative(transcript.start)]
                    )
                    arrows[
                        (
                            loc - spread if strand == "+" else loc + spread,
                            y_loc - exon_width / 5,
                        )
                    ] = Path.MOVETO
                    arrows[(loc, y_loc)] = Path.LINETO
                    arrows[
                        (
                            loc - spread if strand == "+" else loc + spread,
                            y_loc + exon_width / 5,
                        )
                    ] = Path.LINETO

        y_loc += 1
        if plot_domain and obj.domain and transcript.transcript_id in obj.domain.pep:
            y_loc = _plot_domain_pep(
                obj,
                transcript,
                region,
                graph_coords,
                patches,
                y_loc,
                exon_width,
                font_size,
                show_id,
                raster,
                color,
                ax,
            )

    if obj.local_domain:
        y_loc = _plot_local_domain(
            obj,
            region,
            graph_coords,
            patches,
            y_loc,
            exon_width,
            font_size,
            raster,
            color,
            ax,
        )

    ax.add_collection(PatchCollection(patches, match_original=True))
    ax.add_patch(
        PathPatch(
            Path(list(arrows.keys()), list(arrows.values())),
            lw=0.5,
            color=color,
            rasterized=raster,
        )
    )
    ax.set_xlim(min(graph_coords), max(graph_coords))
    ax.set_ylim(-0.5, y_loc + 0.5)


def _plot_domain_pep(
    obj,
    transcript,
    region,
    graph_coords,
    patches,
    y_loc,
    exon_width,
    font_size,
    show_id,
    raster,
    color,
    ax,
):
    """Internal helper for plotting domain protein data."""
    current_domains = obj.domain.pep[transcript.transcript_id]
    for sub_current_domain in current_domains:
        if not (
            sub_current_domain.start <= region.end
            and sub_current_domain.end >= region.start
        ):
            continue
        for sub_exon in sub_current_domain.exons:
            for exon in sub_exon:
                s, e = region.relative(exon.start), region.relative(exon.end)
                if e <= 0 or s > len(region):
                    continue
                s = 0 if s < 0 else s
                e = len(region) - 1 if e > len(region) else e
                patches.append(
                    plt.Rectangle(
                        (graph_coords[s], y_loc - exon_width / 2),
                        width=graph_coords[e] - graph_coords[s],
                        height=exon_width,
                        lw=0.5,
                        zorder=20,
                        rasterized=raster,
                        color=color,
                    )
                )
            intron_relative_s = region.relative(min(map(lambda x_: x_.end, sub_exon)))
            intron_relative_s = intron_relative_s if intron_relative_s >= 0 else 0
            if intron_relative_s > len(region):
                continue
            intron_relative_e = region.relative(max(map(lambda x_: x_.start, sub_exon)))
            intron_relative_e = (
                len(region) - 1
                if intron_relative_e > len(region)
                else intron_relative_e
            )
            if intron_relative_e <= 0:
                continue
            if len(sub_exon) != 1:
                patches.append(
                    plt.Rectangle(
                        (graph_coords[intron_relative_s], y_loc),
                        width=graph_coords[intron_relative_e]
                        - graph_coords[intron_relative_s],
                        height=0,
                        color=color,
                        lw=0.2,
                        rasterized=raster,
                    )
                )
        label_text = (
            f"{sub_current_domain.gene}|{transcript.transcript_id}"
            if show_id
            else f"{sub_current_domain.gene}|{transcript.transcript}"
        )
        ax.text(x=-1, y=y_loc - 0.125, s=label_text, fontsize=font_size / 2, ha="right")
        y_loc += 0.5
    return y_loc


def _plot_local_domain(
    obj, region, graph_coords, patches, y_loc, exon_width, font_size, raster, color, ax
):
    """Internal helper for plotting local domain data."""
    for base_name, current_domain in obj.local_domain.items():
        for sub_current_domain in current_domain:
            if not (
                sub_current_domain.start <= region.end
                and sub_current_domain.end >= region.start
            ):
                continue
            for sub_exon in sub_current_domain.exons:
                for exon in sub_exon:
                    s, e = region.relative(exon.start), region.relative(exon.end)
                    if e <= 0 or s > len(region):
                        continue
                    s = 0 if s < 0 else s
                    e = len(region) - 1 if e > len(region) else e
                    patches.append(
                        plt.Rectangle(
                            (graph_coords[s], y_loc - exon_width / 2),
                            width=graph_coords[e] - graph_coords[s],
                            height=exon_width,
                            lw=0.5,
                            zorder=20,
                            rasterized=raster,
                            color=color,
                        )
                    )
                intron_relative_s = region.relative(
                    min(map(lambda x_: x_.end, sub_exon))
                )
                intron_relative_s = intron_relative_s if intron_relative_s >= 0 else 0
                if intron_relative_s > len(region):
                    continue
                intron_relative_e = region.relative(
                    max(map(lambda x_: x_.start, sub_exon))
                )
                intron_relative_e = (
                    len(region) - 1
                    if intron_relative_e > len(region)
                    else intron_relative_e
                )
                if intron_relative_e <= 0:
                    continue
                if len(sub_exon) != 1:
                    patches.append(
                        plt.Rectangle(
                            (graph_coords[intron_relative_s], y_loc),
                            width=graph_coords[intron_relative_e]
                            - graph_coords[intron_relative_s],
                            height=0,
                            color=color,
                            lw=0.2,
                            rasterized=raster,
                        )
                    )
            ax.text(
                x=-1,
                y=y_loc - 0.125,
                s=f"{sub_current_domain.gene}|{base_name}",
                fontsize=font_size / 2,
                ha="right",
            )
        y_loc += 1
    return y_loc


# ============================================================================
# Density
# ============================================================================


def plot_density(
    ax,
    obj: Optional[File] = None,
    data: Optional[ReadDepth] = None,
    region: Optional[GenomicLoci] = None,
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    color="blue",
    font_size: int = 8,
    show_junction_number: bool = True,
    show_mean_jxn_number: bool = False,
    junction_number_font_size: int = 12,
    n_y_ticks: int = 4,
    distance_between_label_axis: float = 0.1,
    show_y_label: bool = True,
    y_label: str = "",
    theme: str = "ticks_blank",
    max_used_y_val: Optional[float] = None,
    min_used_y_val: Optional[float] = 0,
    raster: bool = False,
    fill_step: str = "post",
    **kwargs,
):
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

    # Junction preparation
    try:
        jxns = data.junctions_dict(show_mean_jxn_number)
    except AttributeError:
        jxns = {}

    jxn_number_log_transformed = False
    if obj and obj.log_trans and obj.log_trans.isdigit() and int(obj.log_trans) > 0:
        jxn_number_log_transformed = True
        y_label += f" (log{obj.log_trans})"
        denominator = np.log(int(obj.log_trans))
        for k, v in jxns.items():
            jxns[k] = np.log1p(v) / denominator

    # Compute data-driven baseline arc heights BEFORE calling _compute_y_limits.
    # These must match the internal computation in _compute_y_limits so that the
    # drawn junction arcs have exactly the same height that was used to expand
    # the y-axis limits.  Using the already-expanded max_used_y_val for arc
    # height causes non-convergent growth and truncated junction arcs.
    if isinstance(data, dict):
        _base_max = max(max(v.plus) if v.plus is not None else 0 for v in data.values())
        _minus_maxes = [
            max(v.minus) if v.minus is not None else 0 for v in data.values()
        ]
        _base_min = -1 * max(_minus_maxes) if _minus_maxes else 0
    else:
        _base_max = max(data.plus) if data.plus is not None else 0
        _base_min = -1 * max(data.minus) if data.minus is not None else 0

    if _base_max % 2 == 1:
        _base_max += 1

    # Arc height reference: use the LARGER of (data baseline, incoming y-limit)
    # so arcs look proportional to the visible axis range (especially in
    # --same-y / --same-y-sc mode).  If we only used the data baseline, panels
    # with small data would have tiny arcs hugging the x-axis when a large
    # global y-limit is shared across panels.
    _arc_ref_max = (
        max(_base_max, max_used_y_val) if max_used_y_val is not None else _base_max
    )
    _arc_ref_min = (
        min(_base_min, min_used_y_val) if min_used_y_val is not None else _base_min
    )

    _top_arc_height = abs(3 * _arc_ref_max / 4)
    _bot_arc_height = (
        abs(3 * _arc_ref_min / 4) if _arc_ref_min != 0 else _top_arc_height
    )

    # Compute y limits using shared logic
    fixed_max_used_y = max_used_y_val is not None

    max_used_y_val, min_used_y_val = _compute_y_limits(
        data=data,
        region=region,
        graph_coords=graph_coords,
        max_used_y_val=max_used_y_val,
        min_used_y_val=min_used_y_val,
        density_by_strand=kwargs.get("density_by_strand", False),
        fill_step=fill_step,
        show_mean_jxn_number=show_mean_jxn_number,
    )
    
    # Draw fill
    x, y1, y2 = [], [], []
    for i in range(len(graph_coords)):
        x.append(graph_coords[i])
        if isinstance(data, dict):
            y1.append(
                max(v.plus[i] if v.plus is not None else 0 for v in data.values())
            )
            y2.append(
                min(-v.minus[i] if v.minus is not None else 0 for v in data.values())
            )
        else:
            y1.append(data.plus[i] if data.plus is not None else 0)
            y2.append(-data.minus[i] if data.minus is not None else 0)

    ax.fill_between(x, y1, y2=y2, color=color, lw=0, step=fill_step, rasterized=raster)

    # Draw junctions
    if jxns:
        jxns_sorted_list = sorted(
            jxns.keys(),
            key=lambda x: (x.end - x.start, x.start, x.end),
            reverse=True,
        )

        max_junction_count = max(jxns.values()) if jxns else 0
        min_junction_count = min(jxns.values()) if jxns else 0
        junction_count_gap = max_junction_count - min_junction_count

        text_objs = []

        for jxn_idx, jxn in enumerate(jxns_sorted_list):
            leftss, rightss = jxn.start, jxn.end
            if (
                (leftss < region.start < region.end < rightss)
                or rightss <= region.start
                or leftss >= region.end
            ):
                logger.debug(
                    f"junction {jxn} of {y_label} is out of plotting region, skip"
                )
                continue
            
            ss1_idx, ss1_modified = get_limited_index(
                leftss - region.start, len(graph_coords)
            )
            ss2_idx, ss2_modified = get_limited_index(
                rightss - region.start, len(graph_coords)
            )
            ss1, ss2 = graph_coords[ss1_idx], graph_coords[ss2_idx]

            if kwargs.get("density_by_strand"):
                jxn_on_top = jxn.strand == "+"
            else:
                jxn_on_top = jxn_idx % 2 == 0
                if abs(min_used_y_val) < max_used_y_val:
                    min_used_y_val = -max_used_y_val

            if fill_step == "post":
                ss1_idx = max(0, ss1_idx - 1)
            elif fill_step == "pre":
                ss2_idx = min(ss2_idx + 1, len(graph_coords) - 1)

            if jxn_on_top:
                left_dens, right_dens = data.curr_max(ss1_idx), data.curr_max(ss2_idx)
                current_height = _top_arc_height
                pts = [
                    left_dens if not ss1_modified else left_dens + current_height,
                    left_dens + current_height,
                    right_dens + current_height,
                    right_dens if not ss2_modified else right_dens + current_height,
                ]

            else:
                left_dens, right_dens = (
                    abs(data.curr_min(ss1_idx)),
                    abs(data.curr_min(ss2_idx)),
                )
                current_height = _bot_arc_height
                pts = [
                    -left_dens if not ss1_modified else -left_dens - current_height,
                    -left_dens - current_height,
                    -right_dens - current_height,
                    -right_dens if not ss2_modified else -right_dens - current_height,
                ]

            line_width = (
                0.5 + 1.5 * (jxns[jxn] - min_junction_count) / junction_count_gap
                if junction_count_gap > 0
                else 0
            )
            pts = [(ss1, pts[0]), (ss1, pts[1]), (ss2, pts[2]), (ss2, pts[3])]

            ax.add_patch(
                PathPatch(
                    Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                    ec=color,
                    lw=line_width + 0.2,
                    fc="none",
                )
            )

            if show_junction_number:
                midpt = cubic_bezier(pts, 0.5)
                val = jxns[jxn]
                val = (
                    round(val, 2)
                    if (jxn_number_log_transformed or show_mean_jxn_number)
                    else int(val)
                )
                t = ax.text(
                    midpt[0],
                    midpt[1],
                    "{0}".format(val),
                    fontsize=junction_number_font_size,
                    ha="center",
                    va="center",
                    backgroundcolor="w",
                )
                t.set_bbox(dict(alpha=0))
                text_objs.append(t)

        if show_junction_number and len(text_objs) > 1:
            adjust_text(
                text_objs,
                ax=ax,
                force_text=(0.5, 0.5),
                expand=(1.2, 1.2),
                arrowprops=None,
                lim=100,
            )

    if obj and obj.title:
        ax.text(
            max(graph_coords) - len(obj.title),
            max_used_y_val,
            obj.title,
            color=color,
            fontsize=font_size,
        )

    if not fixed_max_used_y:
        _, max_used_y_val = ax.get_ylim()
    if min_used_y_val is None:
        min_used_y_val, _ = ax.get_ylim()
    if max_used_y_val == min_used_y_val == 0:
        min_used_y_val, max_used_y_val = ax.get_ylim()

    if (
        not isinstance(data, dict)
        and data.strand_aware
        and kwargs.get("density_by_strand")
    ):
        max_used_y_val = max(abs(min_used_y_val), max_used_y_val)
        min_used_y_val = -max_used_y_val
    elif not kwargs.get("density_by_strand") and not jxns:
        min_used_y_val = 0

    set_y_ticks(
        ax,
        label=y_label,
        theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=max_used_y_val,
        min_used_y_val=min_used_y_val,
        n_y_ticks=n_y_ticks,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        y_axis_skip_zero=False,
    )


# ============================================================================
# Site plot
# ============================================================================


def plot_site_plot(
    ax,
    obj: File,
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    color="blue",
    font_size: int = 8,
    n_y_ticks: int = 3,
    distance_between_label_axis: float = 0.1,
    show_y_label: bool = True,
    y_label: str = "",
    strand_choice: str = None,
    theme="ticks",
    raster: bool = False,
    **kwargs,
):
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
    for label, array_plot in zip(["plus", "minus"], [plus, minus]):
        if strand_choice != "all" and label != strand_choice:
            continue
        array_hist = np.repeat(graph_coords, np.abs(array_plot).astype(np.int32))
        try:
            kde = gaussian_kde(array_hist)
            fit_value = kde.pdf(graph_coords)
        except (ValueError, np.linalg.LinAlgError):
            continue
        fit_value = fit_value / fit_value.max()
        for x, y in zip(graph_coords, array_plot):
            patches.append(
                plt.Rectangle(
                    (x, 0), height=y, width=0.6, color=color, rasterized=raster
                )
            )
        ax.plot(
            graph_coords,
            fit_value * array_plot.max()
            if label == "plus"
            else fit_value * array_plot.min(),
            c=color,
            lw=1,
        )

    ax.add_collection(PatchCollection(patches, match_original=True))
    set_y_ticks(
        ax,
        label=y_label,
        theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=1.1 * max_val,
        min_used_y_val=-1.1 * max_val,
        n_y_ticks=n_y_ticks,
        font_size=font_size,
        show_y_label=False,
        distance_between_label_axis=distance_between_label_axis,
        y_axis_skip_zero=False,
    )


# ============================================================================
# Heatmap
# ============================================================================


def plot_heatmap(
    ax,
    cbar_ax,
    data: Dict[str, ReadDepth],
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    color="viridis",
    font_size: int = 8,
    distance_between_label_axis: float = 0.1,
    show_y_label: bool = True,
    theme: str = "ticks",
    y_label: str = "",
    do_scale: bool = False,
    clustering: bool = False,
    clustering_method: str = "ward",
    distance_metric: str = "euclidean",
    raster: bool = True,
    show_row_names: bool = False,
    vmin=None,
    vmax=None,
    **kwargs,
):
    labels = list(data.keys())
    mtx = np.array([x.wiggle for x in data.values()])

    if clustering and len(mtx) > 1:
        assert clustering_method in CLUSTERING_METHOD, (
            f"clustering_method {clustering_method} is not supported."
        )
        assert distance_metric in DISTANCE_METRIC, (
            f"distance_metric {distance_metric} is not supported"
        )
        order = dendrogram(
            linkage(mtx, method=clustering_method, metric=distance_metric),
            orientation="right",
            no_plot=True,
        )
        mtx = mtx[order["leaves"], :]
        labels = [labels[x] for x in order["leaves"]]

    if do_scale:
        mtx = zscore(mtx, axis=1)

    if not show_row_names:
        labels = False

    if vmin == vmax == 0 and (mtx.max() != 0 or mtx.min() != 0):
        logger.debug("The vmin and vmax == 0, but input matrix not != 0")
        vmin, vmax = None, None

    sns.heatmap(
        mtx,
        ax=ax,
        cmap=color,
        cbar_ax=cbar_ax,
        xticklabels=False,
        yticklabels=labels,
        center=False,
        rasterized=raster,
        vmin=vmin,
        vmax=vmax,
    )

    ax.tick_params(axis="both", which="major", labelsize=font_size, rotation=0)
    cbar_ax.tick_params(labelsize=font_size)

    ymax, ymin = ax.get_ylim()
    set_y_ticks(
        ax,
        label=y_label,
        theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=ymax,
        min_used_y_val=ymin,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        set_label_only=True,
    )


# ============================================================================
# HiC
# ============================================================================


def plot_hic(
    ax,
    cbar_ax,
    obj: List[HiCTrack],
    show_legend: bool = True,
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    color: str = "RdYlBu_r",
    font_size: int = 8,
    distance_between_label_axis: float = 0.1,
    show_y_label: bool = True,
    theme: str = "ticks",
    y_label: str = "",
    n_y_ticks: int = 4,
    raster: bool = True,
    **kwargs,
):
    assert len(obj) == 1, "HiC plot only support one file"
    obj = obj[0]

    if graph_coords is None:
        graph_coords = init_graph_coords(obj.region)
    if not y_label:
        y_label = obj.label

    y_max = obj.matrix.shape[1]
    ax.pcolormesh(
        obj.x - obj.region.start,
        obj.y,
        np.flipud(obj.matrix),
        cmap=color,
        rasterized=raster,
    )

    if obj.tad:
        for tad in obj.tad_list:
            x1 = tad.start - obj.region.start
            x2 = x1 + float(tad.end - tad.start) / 2
            x3 = tad.end - obj.region.start
            y1 = 0
            y2 = x3 - x1
            if y2 > obj.depth:
                continue
            triangle = Polygon(
                [[x1, y1], [x2, y2], [x3, y1]],
                closed=False,
                edgecolor="grey",
                facecolor="None",
            )
            ax.add_artist(triangle)

    ax.set_xlim(0, len(obj.region))
    ax.set_ylim(0, y_max)

    if show_legend:
        cbar = plt.colorbar(
            mappable=mpl.cm.ScalarMappable(
                norm=mpl.colors.Normalize(
                    vmin=0, vmax=np.percentile(obj.matrix.diagonal(1), 80)
                ),
                cmap=color,
            ),
            cax=cbar_ax,
        )
        cbar.ax.tick_params(labelsize=font_size)
        if obj.log_trans:
            legend_ticks = cbar.get_ticks().tolist()
            legend_ticks[0] = f"{legend_ticks[0]}\\n{obj.log_trans}"
            cbar.set_ticklabels(legend_ticks)

    ax.axis("off")
    set_y_ticks(
        ax,
        label=y_label,
        theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=obj.depth,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        n_y_ticks=n_y_ticks,
    )


# ============================================================================
# Line
# ============================================================================


def plot_line(
    ax,
    data: Dict[str, ReadDepth],
    graph_coords: Union[Dict, np.ndarray],
    font_size: int = 8,
    distance_between_label_axis: float = 0.1,
    show_y_label: bool = True,
    y_label: str = "",
    line_attrs: Optional[Dict[str, Dict]] = None,
    theme: str = "ticks_blank",
    n_y_ticks: int = 4,
    show_legend: bool = False,
    legend_position: str = "upper right",
    legend_ncol: int = 0,
    max_used_y_val: Optional[float] = None,
    **kwargs,
):
    max_y_val = 0
    legend = []
    distance_to_zero = plt.linspace(0, len(graph_coords), len(data) + 2)
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
                x[idx],
                np.median(y),
                ylab,
                color=attr.get("color", "black"),
                fontsize=font_size,
                ha="center",
                va="center",
                backgroundcolor="w",
            )
            t.set_bbox(dict(alpha=0))
            legend.append(t)

    if show_legend:
        ax.legend(
            loc=legend_position,
            ncol=int(len(data) / 1.5) if legend_ncol <= 0 else legend_ncol,
            fancybox=False,
            shadow=False,
        )
    else:
        try:
            adjust_text(legend, arrowprops=dict(arrowstyle="-", lw=1))
        except IndexError as err:
            logger.debug(err)

    set_y_ticks(
        ax,
        label=y_label,
        theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=max_y_val if max_used_y_val is None else max_used_y_val,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
        n_y_ticks=n_y_ticks,
        **kwargs,
    )


# ============================================================================
# IGV-like
# ============================================================================


def plot_igv_like(
    ax,
    obj: Dict[str, ReadSegment],
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    y_label: str = "",
    exon_color: Optional[str] = None,
    intron_color: Optional[str] = None,
    feature_color: Optional[str] = None,
    exon_width: float = 0.3,
    font_size: int = 8,
    n_y_ticks: int = 1,
    distance_between_label_axis: float = 0.1,
    show_y_label: bool = True,
    theme: str = "ticks_blank",
    raster: bool = False,
    **kwargs,
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
    if obj.meta is not None:
        for c_ind_list in obj.get_index():
            add_plot = False
            for c_ind in c_ind_list:
                c_data = obj.data[c_ind]
                if c_data.start < region.start or c_data.end > region.end:
                    continue

                for exon in c_data.exons:
                    s, e = region.relative(exon.start), region.relative(exon.end)
                    if e < 0 or s > len(region):
                        continue
                    s = 0 if s < 0 else s
                    e = len(region) - 1 if e >= len(region) else e
                    patches.append(
                        plt.Rectangle(
                            (graph_coords[s], y_loc - exon_width),
                            width=graph_coords[e] - graph_coords[s],
                            height=exon_width * 2,
                            facecolor="k" if not exon_color else exon_color,
                            lw=0.5,
                            zorder=20,
                        )
                    )
                    add_plot = add_plot | True

                for intron in c_data.introns:
                    s, e = region.relative(intron.start), region.relative(intron.end)
                    if e < 0 or s > len(region):
                        continue
                    s = 0 if s < 0 else s
                    e = len(region) - 1 if e >= len(region) else e
                    patches.append(
                        plt.Rectangle(
                            (graph_coords[s], y_loc),
                            width=graph_coords[e] - graph_coords[s],
                            height=0,
                            color="#4d4d4d" if not intron_color else intron_color,
                            lw=0.2,
                            rasterized=raster,
                        )
                    )

                for feature in c_data.features:
                    if feature.start == -1:
                        continue
                    s, e = region.relative(feature.start), region.relative(feature.end)
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
                        patches.append(
                            plt.Rectangle(
                                (graph_coords[s], y_loc - exon_width * width_ratio),
                                width=graph_coords[e] - graph_coords[s],
                                height=exon_width * 2,
                                facecolor="blue"
                                if not feature_color
                                else feature_color,
                                lw=0.2,
                                zorder=20,
                                rasterized=raster,
                            )
                        )
                        scatters_x.append(graph_coords[s])
                        scatters_y.append(y_loc + exon_width * width_ratio)
                    else:
                        patches.append(
                            plt.Rectangle(
                                (graph_coords[s], y_loc - exon_width * width_ratio),
                                width=graph_coords[e] - graph_coords[s],
                                height=exon_width * 2,
                                facecolor="red" if not feature_color else feature_color,
                                lw=0.2,
                                zorder=20,
                                rasterized=raster,
                            )
                        )
            if add_plot:
                y_loc += 1

    ax.add_collection(PatchCollection(patches, match_original=True))
    ax.scatter(
        scatters_x,
        scatters_y,
        rasterized=raster,
        color="b" if not feature_color else feature_color,
        s=0.5,
        linewidths=(0,),
    )

    set_y_ticks(
        ax,
        label=y_label,
        theme=theme,
        graph_coords=graph_coords,
        max_used_y_val=y_loc,
        n_y_ticks=n_y_ticks,
        distance_between_label_axis=distance_between_label_axis,
        font_size=font_size,
        show_y_label=show_y_label,
    )


# ============================================================================
# Links
# ============================================================================


def plot_links(
    ax,
    data: List[Stroke],
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    max_y: int = -10,
    **kwargs,
):
    for stroke in sorted(data):
        leftss, rightss = graph_coords[stroke.start], graph_coords[stroke.end]
        step = (rightss - leftss) / 4
        pts = [
            (leftss, 0),
            (leftss + step, max_y),
            (rightss - step, max_y),
            (rightss, 0),
        ]
        a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
        ax.add_patch(PathPatch(a, fc="none", ec=stroke.color, lw=1))

    ax.set_xlim(min(graph_coords), max(graph_coords))
    ax.set_ylim(max_y, 0)
    ax.axis("off")


# ============================================================================
# Motif
# ============================================================================


def plot_motif(
    ax,
    obj,
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    width: float = 0.8,
    colors=None,
    width_per_character: float = 3.5,
    theme: str = "blank",
    **kwargs,
):
    Theme.set_theme(ax, theme)
    data = obj.data
    if colors is None:
        colors = ["#008000", "#cc0000", "#0000cc", "#ffb300"]
        colors = {x: y for x, y in zip(["A", "T", "C", "G"], colors)}

    region = obj.region
    if graph_coords is None:
        graph_coords = init_graph_coords(region)

    if not data:
        logger.info("there is no any motif information to plot")
        return

    start_base = min(data.keys())
    end_base = max(data.keys())
    ymin, ymax = 0, 0
    xmin = graph_coords[start_base - region.start]
    xmax = graph_coords[end_base - region.start] + (1 + width) / 2

    axin_width = width_per_character * (xmax - xmin) / len(graph_coords)
    bbox_to_left = (xmax * 2 - xmin) / len(graph_coords)
    draw_left = False
    if axin_width + bbox_to_left > 1:
        bbox_to_left = max(xmin / len(graph_coords) - axin_width - 0.1, 0)
        draw_left = True

    axins = ax.inset_axes(
        bounds=[bbox_to_left, 0, axin_width, 1], transform=ax.transAxes
    )

    start_site = start_base
    site = 0
    for idx, vals in data.items():
        site = idx - start_site
        init_height_pos, init_height_neg = 0, 0
        for text, height in vals.items():
            text_shape = make_text_elements(
                text,
                x=site + (1 - width) / 2,
                y=init_height_neg if height < 0 else init_height_pos,
                width=width,
                height=height,
                color=colors.get(text, "blue"),
                edgecolor=colors.get(text, "blue"),
            )
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

    y_height = ymax - ymin
    y_center = (ymax + ymin) / 2
    y_scale = 0.25
    ax.plot(
        [xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymin + y_height * y_scale, ymin + y_height * y_scale, ymin],
        "black",
    )
    ax.set_xlim(min(graph_coords), max(graph_coords))
    ax.set_ylim(ymin, ymax)

    if draw_left:
        pos = len(graph_coords) * ((bbox_to_left + axin_width) * 100 + 5) / 100
        x = xmin
    else:
        pos = len(graph_coords) * bbox_to_left
        x = xmax

    ax.arrow(
        x,
        ymin + y_height * y_scale / 2,
        pos - x,
        y_center - (ymin + y_height * y_scale),
        head_width=1,
        fc="k",
        ec="k",
        rasterized=True,
    )

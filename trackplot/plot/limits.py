#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
Y-axis limit computation — shared between precompute and plot_density.

Dependencies: GenomicLoci, ReadDepth, File, coord.init_graph_coords, utils.get_limited_index
"""

from typing import Dict, Optional, Tuple, Union

import numpy as np

from trackplot.base.GenomicLoci import GenomicLoci
from trackplot.base.ReadDepth import ReadDepth
from trackplot.file.File import File
from trackplot.plot.coord import init_graph_coords
from trackplot.plot.utils import get_limited_index


def _compute_y_limits(
    data: ReadDepth,
    region: GenomicLoci,
    graph_coords: np.ndarray,
    max_used_y_val: Optional[float] = None,
    min_used_y_val: Optional[float] = None,
    density_by_strand: bool = False,
    fill_step: str = "post",
    show_mean_jxn_number: bool = False,
) -> Tuple[float, float]:
    """
    Compute y-axis limits from ReadDepth data, including junction arc extents.

    This is the shared core used by both precompute_y_limits() and plot_density()
    to avoid code duplication.

    Returns (max_used_y_val, min_used_y_val).
    """
    # --- Determine data-driven baseline values for arc height ---
    # These baselines capture the "true" data range before junction expansion.
    # Arc height (current_height) MUST be based on these baselines, not on
    # the running max_used_y_val/min_used_y_val which may already have been
    # expanded by previous junctions.  Using expanded values causes arc
    # heights to grow non-convergently when _compute_y_limits is called
    # again later (e.g. plot_density re-using a precomputed limit).
    if isinstance(data, dict):
        base_max = max(max(v.plus) if v.plus is not None else 0 for v in data.values())
        minus_maxes = [
            max(v.minus) if v.minus is not None else 0 for v in data.values()
        ]
        base_min = -1 * max(minus_maxes) if minus_maxes else 0
    else:
        base_max = max(data.plus) if data.plus is not None else 0
        base_min = -1 * max(data.minus) if data.minus is not None else 0

    # Round base_max up to even for the y-axis limit
    if base_max % 2 == 1:
        base_max += 1

    # --- Set initial y-axis limits ---
    if max_used_y_val is None:
        max_used_y_val = base_max
    if min_used_y_val is None:
        min_used_y_val = base_min

    # --- Arc heights: use the LARGER of (data baseline, incoming y-limit) ---
    # In --same-y mode, max_used_y_val may be the global limit (much larger than
    # this panel's own baseline).  We want the arcs to look proportional to the
    # visible axis range, so we take max(baseline, incoming_limit).
    _arc_ref_max = base_max
    if max_used_y_val is not None:
        _arc_ref_max = max(base_max, max_used_y_val)

    _arc_ref_min = base_min
    if min_used_y_val is not None:
        _arc_ref_min = min(
            base_min, min_used_y_val
        )  # both negative → pick the more negative

    top_arc_height = abs(3 * _arc_ref_max / 4)
    bot_arc_height = abs(3 * _arc_ref_min / 4) if _arc_ref_min != 0 else top_arc_height

    # Get junctions
    try:
        jxns = data.junctions_dict(show_mean_jxn_number)
    except AttributeError:
        jxns = {}

    # Process junctions for y-limit expansion
    if jxns:
        jxns_sorted_list = sorted(
            jxns.keys(),
            key=lambda x: (x.end - x.start, x.start, x.end),
            reverse=True,
        )

        for jxn_idx, jxn in enumerate(jxns_sorted_list):
            leftss, rightss = jxn.start, jxn.end

            if (
                (leftss < region.start < region.end < rightss)
                or rightss <= region.start
                or leftss >= region.end
            ):
                continue

            ss1_idx, ss1_modified = get_limited_index(
                leftss - region.start, len(graph_coords)
            )
            ss2_idx, ss2_modified = get_limited_index(
                rightss - region.start, len(graph_coords)
            )

            if density_by_strand:
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
                left_dens = data.curr_max(ss1_idx)
                right_dens = data.curr_max(ss2_idx)
                current_height = top_arc_height
                pts_y = [
                    left_dens if not ss1_modified else left_dens + current_height,
                    left_dens + current_height,
                    right_dens + current_height,
                    right_dens if not ss2_modified else right_dens + current_height,
                ]
            else:
                left_dens = abs(data.curr_min(ss1_idx))
                right_dens = abs(data.curr_min(ss2_idx))
                current_height = bot_arc_height
                pts_y = [
                    -left_dens if not ss1_modified else -left_dens - current_height,
                    -left_dens - current_height,
                    -right_dens - current_height,
                    -right_dens if not ss2_modified else -right_dens - current_height,
                ]

            max_used_y_val = max(max_used_y_val, max(pts_y))
            min_used_y_val = min(min_used_y_val, min(pts_y))

    # Strand-aware adjustment
    if not isinstance(data, dict) and data.strand_aware:
        if density_by_strand:
            max_used_y_val = max(abs(min_used_y_val), max_used_y_val)
            min_used_y_val = -max_used_y_val
    elif not density_by_strand and not jxns:
        min_used_y_val = 0

    return max_used_y_val, min_used_y_val


def precompute_y_limits(
    obj: Optional[File] = None,
    data: Optional[ReadDepth] = None,
    region: Optional[GenomicLoci] = None,
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    max_used_y_val: Optional[float] = None,
    fill_step: str = "post",
    **kwargs,
):
    """
    Precompute y-axis limits from data + junctions.
    Thin wrapper around _compute_y_limits, used by Plot.plot() for same-y scaling.
    """
    if obj:
        assert obj.region is not None, "please load data first"
        region = obj.region
        if graph_coords is None:
            graph_coords = init_graph_coords(region)
    elif not (region and data):
        raise ValueError("please input obj or region and data")

    if data is None:
        data = obj.data

    return _compute_y_limits(
        data=data,
        region=region,
        graph_coords=graph_coords,
        max_used_y_val=max_used_y_val,
        min_used_y_val=None,
        density_by_strand=kwargs.get("density_by_strand", False),
        fill_step=fill_step,
        show_mean_jxn_number=False,
    )

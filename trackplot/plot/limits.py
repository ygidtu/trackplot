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
    junctions_on_top: bool = False,
) -> Tuple[float, float, float]:
    """
    Compute y-axis limits from ReadDepth data, including junction arc extents.

    This is the shared core used by both precompute_y_limits() and plot_density()
    to avoid code duplication.

    Returns (max_used_y_val, min_used_y_val, base_max).
    base_max is the raw data maximum before junction arc expansion, used as
    arc height reference for --same-y cross-panel consistency.
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

    # --- Track whether limits were explicitly provided (from precompute/user) ---
    # When both are given, they already account for junction expansion, so we
    # must NOT expand them again — doing so causes non-convergent arc growth
    # (the arc height computed from the incoming limit is larger than the one
    # used during precomputation, inflating the limit on every re-calculation).
    _limits_explicitly_set = max_used_y_val is not None and min_used_y_val is not None

    # --- Set initial y-axis limits ---
    if max_used_y_val is None:
        max_used_y_val = base_max
    if min_used_y_val is None:
        min_used_y_val = base_min

    # --- Arc heights: always based on the data baseline ---
    # Previously this used max(base_max, max_used_y_val) which caused
    # arc height inflation when max_used_y_val was already expanded by
    # the precompute step.  The inflated arc peak then exceeded the
    # y-axis limit and overflowed the panel.
    # Use base_max / base_min exclusively — the limit expansion logic
    # below still ensures the axis range grows to accommodate the arcs.
    _arc_ref_max = base_max
    _arc_ref_min = base_min

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
            elif junctions_on_top:
                jxn_on_top = True
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

            # Only expand y-limits for junctions when limits were NOT explicitly
            # provided.  Explicit limits (from precompute / --same-y) already
            # account for junction arcs; re-expanding would cause non-convergent
            # growth because arc heights are recomputed from the incoming limit.
            if not _limits_explicitly_set:
                max_used_y_val = max(max_used_y_val, max(pts_y))
                min_used_y_val = min(min_used_y_val, min(pts_y))

    # Strand-aware adjustment
    if not isinstance(data, dict) and data.strand_aware:
        if density_by_strand:
            max_used_y_val = max(abs(min_used_y_val), max_used_y_val)
            min_used_y_val = -max_used_y_val
    elif not density_by_strand and (junctions_on_top or not jxns):
        min_used_y_val = 0

    # Small expansion so the top / bottom arcs have visual breathing room
    # and do not appear cramped against the axis edge.
    max_used_y_val *= 1.1
    min_used_y_val *= 1.1

    return max_used_y_val, min_used_y_val, base_max


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

    _max, _min, _base = _compute_y_limits(
        data=data,
        region=region,
        graph_coords=graph_coords,
        max_used_y_val=max_used_y_val,
        min_used_y_val=None,
        density_by_strand=kwargs.get("density_by_strand", False),
        fill_step=fill_step,
        show_mean_jxn_number=False,
        junctions_on_top=kwargs.get("junctions_on_top", False),
    )

    return _max, _min, _base

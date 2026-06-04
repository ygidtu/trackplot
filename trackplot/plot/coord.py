#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
Coordinate system transformation, axis helpers, and overlay drawing.

Dependencies: GenomicLoci, Theme, utils
"""

import itertools
import math
from typing import Dict, List, Optional, Union

import numpy as np
from loguru import logger
from matplotlib import pylab

from trackplot.anno.theme import Theme
from trackplot.base.GenomicLoci import GenomicLoci
from trackplot.plot.utils import __merge_exons__


# ============================================================================
# Coordinate system
# ============================================================================


def init_graph_coords(
    region: GenomicLoci,
    exons: Optional[List[List[int]]] = None,
    exon_scale=1,
    intron_scale=0.5,
) -> np.ndarray:
    """
    init the default
    :param region: the plot region
    :param exons: list of start and end sites of exons [[start, end], [start, end]]
    :param exon_scale: the scale of exon, default set to 1
    :param intron_scale: the scale of intron, default set to 0.5 -> the intron showed in final will be half size of real
    """
    graph_coords = np.zeros(len(region), dtype=int)

    if exons:
        if intron_scale <= 1:
            for i in range(0, exons[0][0] - region.start):
                graph_coords[i] = (i - 0) * intron_scale
            exons = __merge_exons__(exons)
            for i in range(0, len(exons)):
                exon = exons[i]
                if i > 0:
                    intron = [exons[i - 1][1], exons[i][0]]
                    for j in range(intron[0], intron[1]):
                        if j >= region.start:
                            graph_coords[j - region.start] = (
                                graph_coords[intron[0] - region.start - 1]
                                + (j - intron[0] + 1) * intron_scale
                            )
                for j in range(exon[0], exon[1] + 1):
                    if j >= region.start:
                        graph_coords[j - region.start] = (
                            graph_coords[exon[0] - region.start - 1]
                            + (j - exon[0] + 1) * exon_scale
                        )
            intron = [exons[-1][-1], region.end]
            for i in range(intron[0], intron[1]):
                if i >= region.start:
                    graph_coords[i - region.start] = (
                        graph_coords[intron[0] - region.start - 1]
                        + (i - intron[0] + 1) * intron_scale
                    )
        else:
            exons = __merge_exons__(exons)
            while exons[0][1] < region.start:
                exons = exons[1:]
            while exons[-1][0] > region.end:
                exons = exons[:-1]
            exons[0][0] = max(exons[0][0], region.start)
            exons[-1][1] = min(exons[-1][1], region.end)

            for i, e in enumerate(exons):
                exons[i][0] -= region.start
                exons[i][1] -= region.start

            steps = [float(exon_scale)] * len(region)
            if exons[0][0] > 0:
                step = intron_scale / exons[0][0]
                for i in range(exons[0][0]):
                    steps[i] = step
            for e in range(1, len(exons)):
                interval = exons[e][0] - exons[e - 1][1] - 2
                step = intron_scale / interval
                for i in range(exons[e - 1][1], exons[e][0]):
                    steps[i] = step

            if last_interval := len(region) - exons[-1][1] - 1:
                step = intron_scale / last_interval
                for i in range(exons[1][1] + 1, len(region)):
                    steps[i] = step
            graph_coords = list(map(int, itertools.accumulate(steps)))
    else:
        for i, j in enumerate(range(region.start, region.end + 1)):
            graph_coords[j - region.start] = i

    if graph_coords[-1] == 0:
        current_max = max(np.where(graph_coords == np.max(graph_coords))[-1])
        for i in np.where(graph_coords == 0)[0]:
            if i > current_max:
                graph_coords[i] = max(graph_coords) + 1
    return graph_coords


# ============================================================================
# Axis helpers
# ============================================================================


def set_x_ticks(
    ax,
    region: GenomicLoci,
    graph_coords: Optional[Union[Dict, np.ndarray]] = None,
    sequence: Optional[Dict[int, str]] = None,
    log_trans: Optional[str] = None,
    nx_ticks: int = 4,
    font_size: int = 6,
    **kwargs,
):
    Theme.set_theme(ax, "blank")
    if graph_coords is None:
        graph_coords = init_graph_coords(region)

    x_label = str(region)

    if log_trans:
        if log_trans in ["2", "10"]:
            log_trans = f"log{log_trans}"
        x_label = f"{x_label}, y axis is {log_trans} transformed"

    ax.hlines(y=0, xmin=0, xmax=max(graph_coords), color="black", lw=1)
    ax.text(
        x=graph_coords[len(graph_coords) // 2],
        y=-2.8,
        s=x_label,
        fontsize=font_size,
        ha="center",
        va="top",
    )

    bk = 1
    if not sequence and nx_ticks > 1:
        bk = max(graph_coords) // (nx_ticks - 1)

    reverse_graph_coords = {}
    for i in range(1, len(graph_coords)):
        for x, y in zip(
            range(graph_coords[i - 1], graph_coords[i]),
            np.linspace(i - 1, i, num=graph_coords[i] - graph_coords[i - 1]),
        ):
            reverse_graph_coords[x] = int(y)

    for i in range(graph_coords[0]):
        reverse_graph_coords[i] = graph_coords[0]

    line_space = {}
    for i in range(nx_ticks - 1):
        i = bk * i
        line_space[i] = reverse_graph_coords[i] + region.start
    line_space[max(graph_coords)] = region.end

    if sequence:
        for i, seq in sequence.items():
            relative_i = graph_coords[i - region.start]
            if relative_i in line_space.keys():
                temp_txs = "{}\\n{}".format(sequence.get(i, ""), line_space[relative_i])
            else:
                temp_txs = "\\n{}".format(sequence.get(i, ""))
            line_space[relative_i] = temp_txs

    for x, s in line_space.items():
        ax.vlines(x=x, ymin=-0.5, ymax=0, color="black", lw=1)
        ax.text(x=x, y=-1, s=s, fontsize=font_size, ha="center", va="top")

    ax.set_ylim(-3.5, 1)
    ax.set_xlim(min(graph_coords), max(graph_coords))


def set_y_ticks(
    ax,
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
    **kwargs,
):
    """
    The y ticks are formatted here
    @2019.03.31 add little check here to make sure the y-axis shows the real value
    """
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

        assign_ticks_y = [
            int(x / (abs(min_used_y_val) + max_used_y_val) * n_y_ticks)
            for x in [abs(min_used_y_val), abs(max_used_y_val)]
        ]
        universal_y_ticks = pylab.linspace(min_used_y_val, 0, assign_ticks_y[0] + 1)
        for i in pylab.linspace(0, max_used_y_val, assign_ticks_y[1] + 1):
            universal_y_ticks = np.append(universal_y_ticks, i)

        universal_y_ticks = np.unique(np.append(universal_y_ticks, 0))
        universal_y_ticks = sorted(universal_y_ticks)

        for lab in universal_y_ticks:
            curr_y_tick_labels.append(f"{int(lab)}")

        ax.set_yticks(universal_y_ticks)
        ax.set_yticklabels(curr_y_tick_labels, fontsize=font_size)
        ax.yaxis.set_ticks_position("left")

    if show_y_label:

        def __dynamic_distance__(
            distance_between_label_axis: float, label: str, scale: int = 100
        ) -> float:
            if distance_between_label_axis != 0:
                return -distance_between_label_axis
            return max(0.01, math.ceil(len(label) / 10) * 10 / scale) * -1

        curr_y_tick_labels = sorted(
            curr_y_tick_labels, key=lambda x: len(x), reverse=True
        )
        ax.text(
            x=__dynamic_distance__(
                distance_between_label_axis,
                curr_y_tick_labels[0] if curr_y_tick_labels else "",
            )
            * max(graph_coords),
            y=(max_used_y_val + min_used_y_val) / 2,
            s=label,
            fontsize=font_size,
            ha="right",
        )

    if max_used_y_val is not None and min_used_y_val is not None:
        ax.set_ylim(ymin=min_used_y_val, ymax=max_used_y_val)


# ============================================================================
# Overlays
# ============================================================================


def set_focus(ax, graph_coords: Union[Dict, np.array], focus: Dict[int, int]):
    for left, right in focus.items():
        try:
            left, right = graph_coords[left], graph_coords[right]
            fill_x = [left, right, right, left]
            y1, y2 = ax.get_ylim()
            fill_y = [y1, y1, y2, y2]
            ax.fill(fill_x, fill_y, alpha=0.1, color="grey")
        except IndexError as err:
            logger.debug("focus region is out of bound: " + str(err))


def set_indicator_lines(
    ax,
    graph_coords: Union[Dict, np.array],
    sites: Dict[int, str],
    min_y_used: Union[int, float] = 0,
    max_y_used: Union[int, float] = None,
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
                lw=0.5,
            )
        except IndexError as err:
            logger.debug("Indicator line is out of bound: " + str(err))

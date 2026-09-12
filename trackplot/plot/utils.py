#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
Pure utility functions — no dependency on other plot/ submodules.
"""

import gzip
from copy import deepcopy
from typing import Dict, List, Set, Tuple

import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D


def load_barcodes(
    barcode: str,
) -> Tuple[Dict[str, Dict[str, Set[str]]], Dict[str, str]]:
    """
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
    """
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
    return (
        p0 * (1 - t) ** 3
        + 3 * t * p1 * (1 - t) ** 2
        + 3 * t**2 * (1 - t) * p2
        + t**3 * p3
    )


def __merge_exons__(exons: List[List[int]]):
    """
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


def make_text_elements(
    text,
    x=0.0,
    y=0.0,
    width=1.0,
    height=1.0,
    color="blue",
    edgecolor="black",
    font=FontProperties(family="monospace"),
):
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

#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
PlotInfo — data container that holds file objects and their metadata for a single track.

Dependencies: File, GenomicLoci, ReadDepth, ReadSegment
"""

from copy import deepcopy
from multiprocessing import Pool
from typing import Dict, List, Optional, Tuple, Union

from trackplot.base.GenomicLoci import GenomicLoci
from trackplot.base.ReadDepth import ReadDepth
from trackplot.file.File import File
from trackplot.file.ReadSegments import ReadSegment


def __load__(args):
    """Helper for multiprocessing data loading."""
    p, args, kwargs = args
    p.load(*args, **kwargs)
    return p


class PlotInfo(object):
    """
    This class is used to collect all the plot information for a single track.
    """

    __slots__ = ["obj", "group", "type", "category", "height_scale"]

    def __init__(self, obj: File, category: str = "", type_: str = "", group: str = ""):
        """
        :param obj: the input file object
        :param category: the input file category
        :param type_: the plot type
        :param group: group name of input file, only used for heatmap/line plots
        """
        self.obj = [obj]
        self.group = group
        self.type = type_
        self.category = [category]
        self.height_scale = None

    def __str__(self) -> str:
        return ";".join([obj.path for obj in self.obj])

    def __hash__(self) -> int:
        if self.type in ["heatmap", "line"]:
            return hash(self.group)
        return hash(
            (
                tuple([hash(x) for x in self.obj]),
                self.group,
                self.type,
                tuple(self.category),
            )
        )

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __ne__(self, other):
        return hash(self) != hash(other)

    def __add__(self, other):
        if self.category == ["heatmap", "line"]:
            self.obj[0] += other.obj
        else:
            self.obj += other.obj
        return self

    @property
    def is_single_cell(self):
        try:
            return self.obj[0].is_single_cell
        except AttributeError:
            return False

    @property
    def data(self) -> Dict[str, Union[ReadDepth, ReadSegment]]:
        data = {}
        for obj in self.obj:
            if isinstance(obj.data, dict):
                data.update(obj.data)
            elif isinstance(obj, ReadSegment):
                data[obj.label] = obj
            else:
                data[obj.label] = obj.data
        return data

    def len(
        self,
        scale: Union[int, float] = 0.25,
        sc_height_ratio: Optional[Dict[str, float]] = None,
        igv_scale: Optional[float] = None,
    ) -> Tuple[int, List]:
        """
        Calculate the number of rows this plot occupies and the height ratios.
        """
        n = 0
        height_ratio = [1]

        if self.is_single_cell:
            if sc_height_ratio is None:
                sc_height_ratio = {"density": 0.2, "heatmap": 0.2}
            height_ratio = [sc_height_ratio.get(self.type, 1)]

        if not self.category:
            pass
        elif self.type == "site-plot" and self.category[0] == "bam":
            n += 2 if not isinstance(self.obj[0], Depth) else 2 * len(self.obj[0])
        elif self.type == "igv":
            read_scale = self.height_scale
            read_scale = igv_scale if read_scale is None else read_scale
            n += self.obj[0].len(scale / 8 if read_scale is None else read_scale)
        elif isinstance(self.obj[0], Depth):
            n += len(self.obj[0])
        else:
            n += 1

        return int(n), [height_ratio[0] for _ in range(n)]

    def add(self, obj: File, category: str = "", type_: str = ""):
        """
        Add a new input file to this group (only for heatmap/line).
        """
        assert type_ in ["heatmap", "line"], "only heatmap/line plot needs add"
        self.obj.append(obj)
        self.category.append(category)
        return self

    def load(self, region: GenomicLoci, n_jobs: int = 0, *args, **kwargs):
        """
        Load data for all file objects in this plot info.
        Supports multiprocessing when n_jobs > 1.
        """
        if n_jobs <= 1:
            for obj in self.obj:
                obj.load(region=region, *args, **kwargs)
        else:
            with Pool(n_jobs) as p:
                kwargs_ = deepcopy(kwargs)
                kwargs_["region"] = region
                self.obj = p.map(__load__, [[o, args, kwargs_] for o in self.obj])

        for obj in self.obj:
            obj.transform()

        if self.type == "density":
            if len(self.obj) > 1:
                data = self.obj[0].data
                for obj in self.obj[1:]:
                    data += obj.data
                for obj in self.obj:
                    obj.data = data

        return self


# Avoid circular import: Depth only needed in len()
from trackplot.file.Depth import Depth  # noqa: E402

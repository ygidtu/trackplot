#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
This file contains the object to handle bam file related issues.

changelog:
    1. add library parameter for determining of read strand at 2022.4.28.

"""
import os
from copy import deepcopy
from typing import Dict, Optional, Set

import numpy as np
from loguru import logger

from trackplot.base.GenomicLoci import GenomicLoci
from trackplot.base.ReadDepth import ReadDepth
from trackplot.base.Readder import Reader
from trackplot.file.File import SingleCell


class ATAC(SingleCell):

    __slots__ = "title", "label", "size_factor", "barcode_groups", "barcode"

    def __init__(self,
                 path: str, barcode_groups: Dict[str, Set[str]], barcode: str, size_factor,
                 label: str = "", title: str = ""):
        u"""
        init this object
        :param label: the left axis label
        :param title: the default title to show in the upper-right of density plot
        :param barcode_groups:
        :param barcode: key of barcode_groups
        :param size_factor
        """
        super().__init__(path, barcode_groups[barcode])
        self.title = title
        self.label = label if label else os.path.basename(path).replace(".bam", "")

        self.size_factor = size_factor
        self.barcode_groups = barcode_groups
        self.barcode = barcode

    @classmethod
    def index(cls, path: str, barcodes: Optional[Dict[str, Set[str]]] = None):
        reverse_barcode_groups = {}  # dict, key => barcode, value => group
        for key, vals in barcodes.items():
            for val in vals:
                reverse_barcode_groups[val.strip()] = key.strip()
        size_factors, sizes = {}, {}
        for values in Reader.read_depth(path=path):
            values = values.split()
            count = int(values[-1])
            barcode = values[-2]
            if barcode not in reverse_barcode_groups.keys():
                continue

            key = reverse_barcode_groups[barcode]

            sizes[key] = sizes.get(key, 0) + 1
            size_factors[key] = size_factors.get(key, 0) + count

        for key, val in size_factors.items():
            size_factors[key] = size_factors[key] * sizes[key]

        del sizes
        median_size_factor = np.median(np.array(list(size_factors.values())))
        return {x: y / median_size_factor for x, y in size_factors.items()}

    @classmethod
    def create(cls,
               path: str,
               label: str = "",
               title: str = "",
               barcode_groups: Optional[Dict[str, Set[str]]] = None,
               barcode: Optional[str] = None,
               size_factors=None):
        u"""
        :param path: the path to bam file
        :param label: the left axis label
        :param title: the default title to show in the upper-right of density plot
        :param barcode_groups:
        :param barcode: key of barcode_groups:
        :param size_factors
        :return:
        """

        return cls(
            path=path,
            label=label,
            title=title,
            barcode=barcode,
            barcode_groups=barcode_groups,
            size_factor=size_factors
        )

    def __hash__(self):
        return hash(self.label)

    def __str__(self) -> str:

        temp = []

        for x in [self.title, self.label, self.path]:
            if x is None or x == "":
                x = "None"
            temp.append(str(x))

        return "\t".join(temp)

    def load(self,
             region: GenomicLoci,
             threshold: int = 0,
             required_strand: Optional[str] = None,
             log_trans: Optional[str] = None,
             **kwargs):
        """
            determine_depth determines the coverage at each base between start_coord and end_coord, inclusive.

            bam_file_path is the path to the bam file used to \
            determine the depth and junctions on chromosome between start_coord and end_coord

        return values:
            depth_vector,
            which is a Numpy array which contains the coverage at each base position between start_coord and end_coord

            spanned_junctions, which is a dictionary containing the junctions supported by reads.
            The keys in spanned_junctions are the
                names of the junctions, with the format chromosome:lowerBasePosition-higherBasePosition
        :param region: GenomicLoci object including the region for calculating coverage
        :param threshold: minimums counts of the given splice junction for visualization
        :param reads1: None -> all reads, True -> only R1 kept; False -> only R2 kept
        :param required_strand: None -> all reads, else reads on specific strand
        :param log_trans: should one of {"10": np.log10, "2": np.log2}
        """
        self.region = region
        self.log_trans = log_trans

        smooth_bin = kwargs.get("smooth_bin", 10)
        half_bin = smooth_bin // 2
        region = deepcopy(region)
        region.start -= half_bin
        region.end += half_bin

        depth_vector = np.zeros(len(region), dtype=int)
        try:
            for _, start, end, barcode, count in Reader.read_depth(path=self.path, region=region):
                # filter reads by 10x barcodes
                start, end, count = int(start), int(end), int(count)
                if not self.empty_barcode():
                    if not self.has_barcode(barcode):
                        continue
                depth_vector[max(start - region.start, 0)] += count
                depth_vector[min(end - region.start, len(depth_vector) - 1)] += count
        except IOError as err:
            logger.error('There is no .bam file at {0}'.format(self.path))
            logger.error(err)
        except ValueError as err:
            logger.error(self.path)
            logger.error(err)

        self.data = ReadDepth(depth_vector)
        self.data.normalize(self.size_factor[self.barcode], format_="atac")

        depth_vector = self.data.wiggle
        for i in range(half_bin, len(self.region)):
            depth_vector[i] = np.mean(depth_vector[i-half_bin:i+half_bin])
        depth_vector = depth_vector[half_bin:len(depth_vector) - half_bin]
        self.data = ReadDepth(depth_vector)
        return self


if __name__ == '__main__':
    pass

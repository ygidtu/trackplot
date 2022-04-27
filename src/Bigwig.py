#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2022.03.22 by ygidtu@gmail.com

This is used to create heatmap images for multiple signal files contains in bigWig file
"""
import os
from typing import List

import numpy as np
import pyBigWig
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import zscore

from src.SpliceRegion import SpliceRegion
from src.logger import logger

__DISTANCE_METRIC__ = [
    "braycurtis", "canberra", "chebyshev",
    "cityblock", "correlation", "cosine",
    "dice", "euclidean", "hamming",
    "jaccard", "jensenshannon", "kulsinski",
    "kulczynski1", "mahalanobis", "matching",
    "minkowski", "rogerstanimoto", "russellrao",
    "seuclidean", "sokalmichener", "sokalsneath",
    "sqeuclidean", "yule"
]

__CLUSTERING_METHOD__ = ["single", "complete", "average", "weighted", "centroid", "median", "ward"]


class Bigwig(object):

    def __init__(
            self,
            files: List[str],
            alias: str = "",
            do_scale: bool = False,
            clustering: bool = False,
            clustering_method: str = "ward",
            distance_metric: str = "euclidean",
            color_map: str = "viridis"
    ):
        u"""

        :param files:
        :param alias:
        :param do_scale: whether to scale the matrix
        :param clustering: whether reorder matrix by clustering
        :param clustering_method: same as  scipy.cluster.hierarchy.linkage
        :param distance_metric: same as scipy.spatial.distance.pdist
        :param color_map: used for seaborn.heatmap, see: https://matplotlib.org/3.5.1/tutorials/colors/colormaps.html
                    'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                    'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
                    'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'
        """
        self.files = [x for x in files if os.path.exists(x)]
        self.alias = alias
        self.data = None
        self.extent = []

        self.do_scale = do_scale
        self.clustering = clustering
        self.clustering_method = clustering_method
        self.distance_metric = distance_metric
        self.color = color_map
        self.raster = False

        if self.clustering:
            if self.clustering_method not in __CLUSTERING_METHOD__:
                raise ValueError(f"{self.clustering_method} is not a supported clustering method")

            if self.distance_metric not in __DISTANCE_METRIC__:
                raise ValueError(f"{self.distance_metric} is not a supported distance metric")

    def prepare(self, region: SpliceRegion, ):
        data = []

        for i, f in enumerate(self.files):
            with pyBigWig.open(f) as r:
                try:
                    data.append(r.values(region.chromosome, region.start, region.end + 1))
                except RuntimeError as e:
                    logger.warning(e)

                    logger.info("may be caused by the mismatch of chromosome")
                    if region.chromosome.startswith("chr"):
                        data.append(r.values(region.chromosome.replace("chr", ""), region.start, region.end + 1))
                    else:
                        data.append(r.values("chr" + region.chromosome, region.start, region.end + 1))

        self.data = np.array(data)

        if self.clustering and self.data.shape[0] > 1:
            data = linkage(self.data, method=self.clustering_method, metric=self.distance_metric)
            order = dendrogram(data, orientation='right')
            self.data = self.data[order["leaves"], :]

        if self.do_scale:
            """
            y = (x – mean) / standard_deviation
            """
            # b = (self.data.transpose() - np.mean(self.data, axis=1)) / np.std(self.data, axis=1)
            # self.data = b.transpose()
            self.data = zscore(self.data, axis=1)
            pass


if __name__ == "__main__":
    pass

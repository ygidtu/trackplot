#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Generate object for processing HiC matrix information

pre-process code was re-wrote based
https://github.com/deeptools/pyGenomeTracks/blob/c42e74e725d22269c33718d9f5df11e0c45c7378/pygenometracks/tracks/HiCMatrixTrack.py#L13
"""

import itertools
from typing import Optional

import numpy as np
from loguru import logger
from scipy import sparse

from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.Readder import Reader


class HiCTrack:
    def __init__(self,
                 path: str,
                 label: str = "",
                 depth: int = 30000,
                 trans: Optional[str] = None,
                 matrix: Optional[np.ndarray] = None,
                 x_coord: Optional[np.ndarray] = None,
                 y_coord: Optional[np.ndarray] = None,
                 region: Optional[GenomicLoci] = None,
                 is_single_cell: bool = False
                 ):
        self.path = path
        self.matrix = matrix
        self.x = x_coord
        self.y = y_coord
        self.depth = depth
        self.trans = trans
        self.label = label
        self.region = region
        self.is_single_cell = is_single_cell

    @classmethod
    def create(cls,
               path: str,
               label: str,
               depth: int,
               trans: Optional[str] = False
               ):
        """
        Create a HiCTrack object for fetching interaction matrix
        :param path: the HiC file which could be one of [h5, cool / mcool / scool, hicpro, homer]
        :param label: the label of the given HiC data
        :param depth: the depth of the given HiC data, a bigger depth means big y-axis
        :param trans: log1p, log2 or log10 transform
        :return:
        """
        return cls(
            path=path,
            label=label,
            depth=depth,
            trans=trans
        )

    def load(self,
             region: GenomicLoci,
             **kwargs
             ):
        """
        Load data from the given region
        :param region: the GenomicLoci object of given region
        :return:
        """
        hic = Reader.read_hic(path=self.path, region=region)
        chromosome_id = region.chromosome
        region_start = region.start
        region_end = region.end
        try:
            chr_start_id, chr_end_id = hic.getChrBinRange(chromosome_id)
        except:
            logger.info("may be caused by the mismatch of chromosome")
            if chromosome_id.startswith("chr"):
                chromosome_id = region.chromosome.replace("chr", "") if region.chromosome != "chrM" else "MT"
            else:
                chromosome_id = "chr" + region.chromosome if region.chromosome != "MT" else "chrM"

            chr_start_id, chr_end_id = hic.getChrBinRange(chromosome_id)

        chr_start = hic.cut_intervals[chr_start_id][1]
        chr_end = hic.cut_intervals[chr_end_id - 1][2]

        start_bp = max(chr_start, region_start - self.depth)
        end_bp = min(chr_end, region_end + self.depth)

        idx = [idx for idx, x in enumerate(hic.cut_intervals)
               if x[0] == chromosome_id and x[1] >= start_bp and x[2] <= end_bp]

        start_pos = [x[1] for i, x in enumerate(hic.cut_intervals) if i in idx]
        start_pos = tuple(list(start_pos) + [hic.cut_intervals[idx[-1]][2]])

        matrix = hic.matrix[idx, :][:, idx]
        region_len = region_end - region_start

        current_depth = min(self.depth, int(region_len * 1.25))
        depth_in_bins = max(1, int(1.5 * region_len / hic.getBinSize()))

        if current_depth < self.depth:
            logger.warning(f"The depth was set to {self.depth} which is more than 125% "
                           f"of the region plotted. The depth will be set "
                           f"to {current_depth}.\n")
            # remove from matrix all data points that are not visible.
            matrix = matrix - sparse.triu(matrix, k=depth_in_bins, format='csr')

        matrix = np.asarray(matrix.todense().astype(float))

        n = matrix.shape[0]
        t = np.array([[1, 0.5], [-1, 0.5]])
        matrix_tmp = np.dot(np.array([(i[1], i[0])
                                      for i in itertools.product(start_pos[::-1],
                                                                 start_pos)]), t)

        self.x = matrix_tmp[:, 1].reshape(n + 1, n + 1)
        self.y = matrix_tmp[:, 0].reshape(n + 1, n + 1)
        self.matrix = HiCTrack.mat_trans(matrix, method=self.trans)
        self.region = region

    @staticmethod
    def mat_trans(matrix: np.ndarray, method: Optional[str] = None):
        if not method:
            return matrix
        assert method in {"log1p", "log2", "log10"}, f"couldn't recognize current transform method {method}"
        if method == "log1p":
            return np.log1p(matrix)
        elif method == "log2":
            return np.log2(matrix + 1)
        elif method == "log10":
            return np.log10(matrix + 1)


if __name__ == '__main__':
    pass

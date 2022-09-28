#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
This file contains the object to handle bam file related issues.

changelog:
    1. add library parameter for determining of read strand at 2022.4.28.

"""
import gzip
import os
from typing import Optional, Set

import numpy as np
from loguru import logger

from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.ReadDepth import ReadDepth
from sashimi.base.Readder import Reader
from sashimi.file.File import SingleCell


class ATAC(SingleCell):
    def __init__(self, path: str, label: str = "", title: str = "", barcodes: Optional[Set[str]] = None):
        u"""
        init this object
        :param label: the left axis label
        :param title: the default title to show in the upper-right of density plot
        :param barcodes: the path to barcodes,
                default: ../filtered_feature_bc_matrix/barcodes.tsv.gz of bam file according to 10X Genomics
        """
        super().__init__(path, barcodes)
        self.title = title
        self.label = label if label else os.path.basename(path).replace(".bam", "")

    @classmethod
    def create(cls,
               path: str,
               label: str = "",
               title: str = "",
               barcodes: Optional[Set[str]] = None
               ):
        u"""

        :param path: the path to bam file
        :param label: the left axis label
        :param title: the default title to show in the upper-right of density plot
        :param barcodes: the path to barcodes,
                default: ../filtered_feature_bc_matrix/barcodes.tsv.gz of bam file according to 10X Genomics
        :return:
        """
        barcode = barcodes
        path = os.path.abspath(path)
        if not barcodes:
            barcode = set()
            barcodes = os.path.join(os.path.dirname(path), "filtered_feature_bc_matrix/barcodes.tsv.gz")

            if os.path.exists(barcodes):
                with gzip.open(barcodes, "rt") as r:
                    for line in r:
                        barcode.add(line.strip())
        return cls(
            path=path,
            label=label,
            title=title,
            barcodes=barcode
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
             reads1: Optional[bool] = None,
             required_strand: Optional[str] = None,
             log_trans: Optional[str] = None,
             **kwargs
             ):
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

        depth_vector = np.zeros(len(region), dtype=int)

        try:
            for _, start, end, barcode, count in Reader.read_depth(path=self.path, region=region):
                # filter reads by 10x barcodes
                if not self.empty_barcode():
                    if not self.has_barcode(barcode):
                        continue
        except IOError as err:
            logger.error('There is no .bam file at {0}'.format(self.path))
            logger.error(err)
        except ValueError as err:
            logger.error(self.path)
            logger.error(err)

        self.data = ReadDepth(depth_vector)
        return self


if __name__ == '__main__':
    bam = ATAC.create("../../example/bams/sc.bam")
    bam.load(GenomicLoci("chr1", 1270656, 1284730, "+"), 10)

    print(str(bam))
    pass

#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2022.03.22 by ygidtu@gmail.com

This is used to create heatmap images for multiple signal files contains in bigWig file
"""
import os

import numpy as np

from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.ReadDepth import ReadDepth
from sashimi.base.Readder import Reader
from sashimi.file.File import File


class Bigwig(File):

    def __init__(self, path: str, label: str = "", title: str = ""):
        u"""
        :param path: the path to bam file
        :param label: the left axis label
        """
        super().__init__(path)
        self.label = label
        self.title = title

    @classmethod
    def create(cls, path: str, label: str = "", title: str = ""):
        assert os.path.exists(path), f"{path} is not exists."
        if not label:
            label = os.path.basename(path)
        return cls(path=path, label=label, title=title)

    def load(self, region: GenomicLoci, **kwargs):
        self.region = region

        vals = np.nan_to_num(
            Reader.read_bigwig(self.path, region),
            copy=True, nan=0.0, posinf=None, neginf=None
        )
        plus, minus = np.zeros(len(vals)), np.zeros(len(vals))
        for idx, val in enumerate(vals):
            if val > 0:
                plus[idx] = val
            else:
                minus[idx] = abs(val)

        self.data = ReadDepth(plus, minus=minus)


if __name__ == "__main__":
    pass

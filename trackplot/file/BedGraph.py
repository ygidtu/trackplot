#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2022.03.22 by ygidtu@gmail.com

This is used to create heatmap images for multiple signal files contains in bigWig file
"""
import os

import numpy as np

from trackplot.base.GenomicLoci import GenomicLoci
from trackplot.base.ReadDepth import ReadDepth
from trackplot.base.Readder import Reader
from trackplot.file.File import File


class Bedgraph(File):

    __slots__ = "label", "title"

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

        depth = np.zeros(len(self.region))

        for record in Reader.read_depth(self.path, self.region):
            if len(record) >= 4:
                _, start, end, cov = record[:4]
                start, end, cov = int(start), int(end), float(cov)
                depth[(start - self.region.start):(end - self.region.start)] += cov
        self.data = ReadDepth(depth)


if __name__ == "__main__":
    pass

#!/usr/bin/pyton3
# -*- coding:utf-8 -*-
u"""

"""
import os

from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.Readder import Reader


class Motif(object):
    def __init__(self,
                 path: str,
                 region: GenomicLoci,
                 ):
        u"""
        the plot region
        """
        self.path = path
        self.region = region
        self.label = ""
        self.data = {}  # List[List[base, score]]

    @classmethod
    def create(cls, path: str, region: GenomicLoci):
        assert os.path.exists(path), f"{path} is not exists."
        return cls(path=path, region=region)

    def load(self, region: GenomicLoci, **kwargs):
        data = {}
        keys = ["A", "T", "C", "G"]
        for record in Reader.read_depth(self.path, region):
            start = record[1]
            data[int(start)] = {x: float(y) for x, y in zip(keys, record[3:7])}
        self.data = data


if __name__ == '__main__':
    pass

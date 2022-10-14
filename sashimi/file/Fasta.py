#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
u"""

"""
import os

import pysam
from loguru import logger

from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.file.File import File


class Fasta(File):

    def __init__(self, path: str):
        super().__init__(path)
        assert os.path.exists(path), f"{path} not exists"

    @classmethod
    def create(cls, path: str):
        if not os.path.exists(path + ".fai"):
            logger.warning(f"{path}.fai not exists, try to create it")

            pysam.faidx(path)

        return cls(path)

    def load(self, region: GenomicLoci, **kwargs):
        u"""
        load sequence from a fasta file
        :param region:
        :return:
        """

        if not os.path.exists(self.path + ".fai"):
            logger.warning(f"{self.path}.fai not exists, try to create it")
            pysam.faidx(self.path)

        self.data = {}
        with pysam.FastaFile(self.path) as r:
            for i, j in zip(
                    range(self.region.start, self.region.end + 1),
                    r.fetch(region.chromosome, region.start, region.end + 1)
            ):
                self.data[i] = j
        return self


if __name__ == '__main__':
    pass

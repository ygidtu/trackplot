#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
The parent object of input files
"""


class File(object):
    def __init__(self, path: str):
        self.path = path
        self.data = None
        self.label = ""
        self.region = None
        self.log_trans = "0"
        self.title = ""

    @property
    def chrom(self) -> str:
        if self.region:
            return self.region.chrom
        return ""

    @property
    def start(self) -> int:
        if self.region:
            return self.region.start
        return 0

    @property
    def end(self) -> int:
        if self.region:
            return self.region.end
        return 0

if __name__ == "__main__":
    pass

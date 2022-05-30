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

    def load(self, *args, **kwargs):
        return None

    def __hash__(self) -> int:
        return hash((self.path, self.label, self.title))

    def __eq__(self, other):
        return self.path == other.path and self.label == other.label and self.title == other.title


if __name__ == "__main__":
    pass

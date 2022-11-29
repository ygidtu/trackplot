#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
The parent object of input files
"""
from copy import deepcopy
from typing import Optional, Set, List, Dict


class File(object):
    def __init__(self, path: str):
        self.path = path
        self.data = None
        self.label = ""
        self.region = None
        self.log_trans = "0"
        self.title = ""
        self.is_single_cell = False

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

    def len(self, scale=1) -> int:
        return len(self.data) / scale if self.data else 0

    def __hash__(self) -> int:
        return hash((self.path, self.label, self.title))

    def __eq__(self, other):
        return self.path == other.path and self.label == other.label and self.title == other.title

    def transform(self):
        if self.data is not None:
            self.data.transform(self.log_trans)


def __set_barcodes__(barcodes: Optional[List[str]]) -> Dict:
    u"""
    separate barcodes by its first character to reduce set size
    :params barcodes: list or set of barcodes
    """
    res = {}

    if barcodes is not None:
        for b in barcodes:
            if b:
                f = b[:min(3, len(b))]

                if f not in res.keys():
                    res[f] = set()

                res[f].add(b)

    return res


class SingleCell(File):
    def __init__(self, path: str, barcodes: Optional[Set[str]] = None, barcode_tag: str = "CB", umi_tag: str = "UB"):

        super().__init__(path)
        self.barcodes = __set_barcodes__(barcodes)
        self.barcode_tag = barcode_tag
        self.umi_tag = umi_tag
        self.is_single_cell = not self.empty_barcode()

    def has_barcode(self, barcode: str) -> bool:
        u"""
        check whether contains barcodes
        :param barcode: barcode string
        """
        if barcode:
            f = barcode[:min(3, len(barcode))]

            temp = self.barcodes.get(f, set())

            return barcode in temp
        return False

    def empty_barcode(self) -> bool:
        u"""
        check whether this bam do not contain any barcodes

        """
        count = 0

        for i in self.barcodes.values():
            count += len(i)

            if count > 0:
                return False

        return True

    def __add__(self, other):
        self.path += other.path

        for i, j in other.barcodes.items():
            if i not in self.barcodes.keys():
                self.barcodes[i] = j
            else:
                self.barcodes[i] |= j

        return self

    def __eq__(self, other) -> bool:
        return self.__hash__() == other.__hash__()

    def copy(self):
        return deepcopy(self)


if __name__ == "__main__":
    pass

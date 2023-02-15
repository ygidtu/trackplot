#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06

Changelog:
    1. move several attributes and functions to corresponding file objects, turn this into pure data class
    2. add transform to log2, log10 or zscore transform the data while plotting
"""
from typing import Optional

import numpy as np
from scipy.stats import zscore


class ReadDepth(object):
    u"""
    Migrated from SplicePlot ReadDepth class

    add a parent class to handle all the position comparison
    """

    __slots__ = [
        "junctions_dict", "minus", "plus",
        "strand_aware", "site_plus", "site_minus",
    ]

    def __init__(self,
                 wiggle: np.array,
                 site_plus: Optional[np.array] = None,
                 site_minus: Optional[np.array] = None,
                 minus: Optional[np.array] = None,
                 junction_dict: Optional[np.array] = None,
                 strand_aware: bool = False):
        u"""
        init this class
        :param wiggle: a numpy.ndarray object represented the whole read coverage,
                       should be summation of plus and minus or plus
        :param site_plus: a numpy.ndarray object represented the forward site coverage
        :param site_minus: a numpy.ndarray object represented the reverse site coverage
        :param minus: a numpy.ndarray object represented the reverse strand read coverage
        :param strand_aware: strand specific depth
        """
        self.plus = wiggle
        self.strand_aware = strand_aware
        self.minus = abs(minus) if minus is not None else minus
        self.junctions_dict = junction_dict
        self.site_plus = site_plus
        self.site_minus = site_minus * -1 if site_minus is not None else site_minus

    @property
    def wiggle(self) -> np.array:
        if (self.plus is None or not self.plus.any()) and self.minus is not None:
            return self.minus

        if self.plus is not None and self.minus is not None:
            return self.plus + self.minus

        return self.plus

    @property
    def max(self) -> float:
        return max(self.wiggle, default=0)

    def __add__(self, other):
        """
            __add__ allows two ReadDepth objects to be added together using the + symbol

            Both self and other must have the same low and high attributes

            return value:
                A new ReadDepth object containing the sum of the two original ReadDepth objects
        """

        if self.wiggle is not None and other.wiggle is not None:
            if len(self.wiggle) == len(other.wiggle):
                junc = {}
                junc.update(self.junctions_dict)
                junc.update(other.junctions_dict)

                minus = None
                if self.minus is not None and other.minus is not None:
                    minus = self.minus + other.minus
                elif self.minus is None and other.minus is not None:
                    minus = other.minus
                elif self.minus is not None and other.minus is None:
                    minus = self.minus

                return ReadDepth(
                    self.plus + other.plus, minus=minus,
                    junction_dict=junc,
                )
        elif self.wiggle is None:
            return other
        else:
            return self

    def curr_height(self, pos: int) -> float:
        if self.minus is None:
            return self.plus[pos]
        return self.plus[pos] + self.minus[pos]

    def curr_max(self, pos: int) -> float:
        return self.plus[pos]

    def curr_min(self, pos: int) -> float:
        return self.minus[pos] if self.minus is not None else 0

    def add_customized_junctions(self, other):
        u"""
        Add customized junctions to plot
        :param other:
        :return:
        """

        for key, value in other.junctions_dict.items():
            if key in self.junctions_dict.keys():
                self.junctions_dict[key] = value + other.junctions_dict[key]
            else:
                self.junctions_dict[key] = value

        return self.junctions_dict

    def transform(self, log_trans: str):
        funcs = {"10": np.log10, "2": np.log2, "zscore": zscore, "e": np.log}

        if log_trans in funcs.keys():
            if self.plus is not None:
                self.plus = funcs[log_trans](self.plus + 1)

            if self.minus is not None:
                self.minus = funcs[log_trans](self.minus + 1)

    def normalize(self, size_factor: float, format_: str = "normal", read_length: float = 1):
        u"""
        Convert reads counts to cpm, fpkm or just scale with scale_factor

        Inspired by `rpkm_per_region` from [MISO](https://github.com/yarden/MISO/blob/b71402188000465e3430736a11ea118fd5639a4a/misopy/sam_rpkm.py#L51)
        """

        if format_ == "rpkm":
            # for rpkm the size_factor is total reads
            self.plus = np.divide(self.plus, np.multiply((len(self.plus) - read_length + 1) * 1e3, size_factor / 1e6))
            if self.minus is not None:
                self.minus = np.divide(self.minus,
                                       np.multiply((len(self.minus) - read_length + 1) * 1e3, size_factor / 1e6))
        elif format_ == "cpm":
            # for cpm the size_factor is total reads
            self.plus = np.divide(self.plus, np.multiply(size_factor, 1e6))
            if self.minus is not None:
                self.minus = np.divide(self.minus, np.multiply(size_factor, 1e6))
        else:
            self.plus = np.divide(self.plus, size_factor)  # * 100
            if self.minus is not None:
                self.minus = np.divide(self.minus, size_factor)


if __name__ == '__main__':
    pass

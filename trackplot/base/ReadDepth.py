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
        "_junction_dict_plus_", "_junction_dict_minus_",
        "_minus_", "_plus_", "_number_of_merged_",
        "strand_aware", "site_plus", "site_minus",
    ]

    def __init__(self,
                 wiggle: np.array,
                 site_plus: Optional[np.array] = None,
                 site_minus: Optional[np.array] = None,
                 minus: Optional[np.array] = None,
                 junction_dict_plus: Optional[np.array] = None,
                 junction_dict_minus: Optional[np.array] = None,
                 strand_aware: bool = False):
        u"""
        init this class
        :param wiggle: a numpy.ndarray object represented the whole read coverage,
                       should be summation of plus and minus or plus
        :param site_plus: a numpy.ndarray object represented the forward site coverage
        :param site_minus: a numpy.ndarray object represented the reverse site coverage
        :param minus: a numpy.ndarray object represented the reverse strand read coverage
        :param strand_aware: strand specific depth
        :param junction_dict_plus: these splice junction from plus strand
        :param junction_dict_minus: these splice junction from minus strand
        """
        self._plus_ = wiggle
        self.strand_aware = strand_aware
        self._minus_ = abs(minus) if minus is not None else minus
        self._junction_dict_plus_ = junction_dict_plus
        self._junction_dict_minus_ = junction_dict_minus
        self.site_plus = site_plus
        self.site_minus = site_minus * -1 if site_minus is not None else site_minus

        self._number_of_merged_ = 1

    @property
    def plus(self) -> Optional[np.array]:
        if self._plus_ is not None and self._number_of_merged_ > 0:
            return self._plus_ / self._number_of_merged_
        return self._plus_

    @property
    def minus(self) -> Optional[np.array]:
        if self._minus_ is not None and self._number_of_merged_ > 0:
            return self._minus_ / self._number_of_merged_
        return self._minus_

    @property
    def wiggle(self) -> np.array:
        if (self._plus_ is None or not self._plus_.any()) and self._minus_ is not None:
            return self.minus

        if self._plus_ is not None and self._minus_ is not None:
            return self.plus + self.minus

        return self.plus

    @property
    def junctions_plus(self) -> dict:
        if self._number_of_merged_ > 1:
            return {k: v / self._number_of_merged_ for k, v in self._junction_dict_plus_.items()}
        return self._junction_dict_plus_

    @property
    def junctions_minus(self) -> dict:
        if self._number_of_merged_ > 1:
            return {k: v / self._number_of_merged_ for k, v in self._junction_dict_minus_.items()}
        return self._junction_dict_minus_

    @property
    def junctions_dict(self) -> dict:
        res = {}
        if self._junction_dict_plus_:
            res.update(self.junctions_plus)

        if self._junction_dict_minus_:
            res.update(self.junctions_minus)
        return res

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
                junc_plus, junc_minus = {}, {}

                for i in [self._junction_dict_plus_, other._junction_dict_plus_]:
                    if i:
                        junc_plus.update(i)
                for i in [self._junction_dict_minus_, other._junction_dict_minus_]:
                    if i:
                        junc_minus.update(i)

                minus = None
                if self._minus_ is not None and other._minus_ is not None:
                    minus = self._minus_ + other._minus_
                elif self._minus_ is None and other._minus_ is not None:
                    minus = other.minus
                elif self._minus_ is not None and other._minus_ is None:
                    minus = self._minus_

                merged = ReadDepth(
                    self._plus_ + other._plus_, minus=minus,
                    junction_dict_plus=junc_plus,
                    junction_dict_minus=junc_minus
                )
                merged._number_of_merged_ = self._number_of_merged_ + other._number_of_merged_
                return merged
            else:
                raise ValueError(f"ReadDepth objects are not equal length: {len(self.wiggle)} != {len(other.wiggle)}")
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

        for k, v in other._junction_dict_plus_:
            self._junction_dict_plus_[k] = v + self._junction_dict_plus_.get(k, 0)

        for k, v in other._junction_dict_minus_:
            self._junction_dict_minus_[k] = v + self._junction_dict_minus_.get(k, 0)

        return self.junctions_dict

    def transform(self, log_trans: str):
        funcs = {"10": np.log10, "2": np.log2, "zscore": zscore, "e": np.log}

        if log_trans in funcs.keys():
            if self._plus_ is not None:
                self._plus_ = funcs[log_trans](self._plus_ + 1)

            if self.minus is not None:
                self._minus_ = funcs[log_trans](self._minus_ + 1)

    def normalize(self, size_factor: float, format_: str = "normal", read_length: float = 0):
        u"""
        Convert reads counts to cpm, fpkm or just scale with scale_factor

        Inspired by `rpkm_per_region` from
        [MISO](https://github.com/yarden/MISO/blob/b71402188000465e3430736a11ea118fd5639a4a/misopy/sam_rpkm.py#L51)
        """

        if format_ == "rpkm" and read_length > 0:
            # for rpkm the size_factor is total reads
            self._plus_ = np.divide(
                self._plus_,
                np.multiply(
                    (np.sum(self._plus_ != 0) - read_length + 1) / 1e3,
                    size_factor / 1e6
                )
            )
            if self._minus_ is not None:
                self._minus_ = np.divide(
                    self._minus_,
                    np.multiply(
                        (np.sum(self._minus_ != 0) - read_length + 1) / 1e3,
                        size_factor / 1e6
                    )
                )
            elif format_ == "cpm" and read_length > 0:
                # for cpm the size_factor is total reads
                self._plus_ = np.divide(self._plus_, np.divide(size_factor, 1e6))
                if self._minus_ is not None:
                    self._minus_ = np.divide(self._minus_, np.divide(size_factor, 1e6))
        elif format_ == "cpm" and read_length > 0:
            # for cpm the size_factor is total reads
            self._plus_ = np.divide(self._plus_, np.divide(size_factor, 1e6))
            if self._minus_ is not None:
                self._minus_ = np.divide(self._minus_, np.divide(size_factor, 1e6))
        elif size_factor is not None and size_factor > 0 and format_ == "atac":
            self._plus_ = np.divide(self._plus_, size_factor)  # * 100
            if self._minus_ is not None:
                self._minus_ = np.divide(self._minus_, size_factor)


if __name__ == '__main__':
    pass

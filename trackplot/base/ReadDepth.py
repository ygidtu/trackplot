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
        "__junction_dict__plus___", "__junction_dict__minus__",
        "__minus__", "__plus__", "__number_of_merged__",
        "strand_aware", "site_plus", "site_minus",
    ]

    def __init__(self,
                 wiggle: np.ndarray,
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
        self.__plus__ = wiggle
        self.strand_aware = strand_aware
        self.__minus__ = abs(minus) if minus is not None else minus
        self.__junction_dict__plus___ = junction_dict_plus
        self.__junction_dict__minus__ = junction_dict_minus
        self.site_plus = site_plus
        self.site_minus = site_minus * -1 if site_minus is not None else site_minus

        self.__number_of_merged__ = 1

    @property
    def plus(self) -> Optional[np.array]:
        if self.__plus__ is not None and self.__number_of_merged__ > 0:
            return self.__plus__ / self.__number_of_merged__
        return self.__plus__

    @property
    def minus(self) -> Optional[np.array]:
        if self.__minus__ is not None and self.__number_of_merged__ > 0:
            return self.__minus__ / self.__number_of_merged__
        return self.__minus__

    @property
    def wiggle(self) -> np.ndarray:
        if (self.__plus__ is None or not self.__plus__.any()) and self.__minus__ is not None:
            return self.minus

        if self.__plus__ is not None and self.__minus__ is not None:
            return self.plus + self.minus

        return self.plus

    @property
    def mean_junctions_plus(self) -> dict:
        if self.__number_of_merged__ > 1:
            return {k: v / self.__number_of_merged__ for k, v in self.__junction_dict__plus___.items()}
        return self.__junction_dict__plus___

    @property
    def mean_junctions_minus(self) -> dict:
        if self.__number_of_merged__ > 1:
            return {k: v / self.__number_of_merged__ for k, v in self.__junction_dict__minus__.items()}
        return self.__junction_dict__minus__

    def junctions_dict(self, show_mean_jxn_number: bool = False) -> dict:
        res = {}

        if show_mean_jxn_number:
            if self.__junction_dict__plus___:
                res.update(self.mean_junctions_plus)

            if self.__junction_dict__minus__:
                res.update(self.mean_junctions_minus)
        else:
            if self.__junction_dict__plus___:
                res.update(self.__junction_dict__plus___)

            if self.__junction_dict__minus__:
                res.update(self.__junction_dict__minus__)
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

                for i in [self.__junction_dict__plus___, other.__junction_dict__plus___]:
                    if i:
                        junc_plus.update(i)
                for i in [self.__junction_dict__minus__, other.__junction_dict__minus__]:
                    if i:
                        junc_minus.update(i)

                minus = None
                if self.__minus__ is not None and other.__minus__ is not None:
                    minus = self.__minus__ + other.__minus__
                elif self.__minus__ is None and other.__minus__ is not None:
                    minus = other.minus
                elif self.__minus__ is not None and other.__minus__ is None:
                    minus = self.__minus__

                merged = ReadDepth(
                    self.__plus__ + other.__plus__, minus=minus,
                    junction_dict_plus=junc_plus,
                    junction_dict_minus=junc_minus
                )
                merged.__number_of_merged__ = self.__number_of_merged__ + other.__number_of_merged__
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

        for k, v in other.__junction_dict__plus___:
            self.__junction_dict__plus___[k] = v + self.__junction_dict__plus___.get(k, 0)

        for k, v in other.__junction_dict__minus__:
            self.__junction_dict__minus__[k] = v + self.__junction_dict__minus__.get(k, 0)

        return self.junctions_dict

    def transform(self, log_trans: str):
        funcs = {"10": np.log10, "2": np.log2, "zscore": zscore, "e": np.log}

        if log_trans in funcs.keys():
            if self.__plus__ is not None:
                self.__plus__ = funcs[log_trans](self.__plus__ + 1)

            if self.minus is not None:
                self.__minus__ = funcs[log_trans](self.__minus__ + 1)

    def normalize(self, size_factor: float, format_: str = "normal", read_length: float = 0):
        u"""
        Convert reads counts to cpm, fpkm or just scale with scale_factor

        Inspired by `rpkm_per_region` from
        [MISO](https://github.com/yarden/MISO/blob/b71402188000465e3430736a11ea118fd5639a4a/misopy/sam_rpkm.py#L51)
        """

        if format_ == "rpkm" and read_length > 0:
            # for rpkm the size_factor is total reads
            self.__plus__ = np.divide(
                self.__plus__,
                np.multiply(
                    (np.sum(self.__plus__ != 0) - read_length + 1) / 1e3,
                    size_factor / 1e6
                )
            )
            if self.__minus__ is not None:
                self.__minus__ = np.divide(
                    self.__minus__,
                    np.multiply(
                        (np.sum(self.__minus__ != 0) - read_length + 1) / 1e3,
                        size_factor / 1e6
                    )
                )
            elif format_ == "cpm" and read_length > 0:
                # for cpm the size_factor is total reads
                self.__plus__ = np.divide(self.__plus__, np.divide(size_factor, 1e6))
                if self.__minus__ is not None:
                    self.__minus__ = np.divide(self.__minus__, np.divide(size_factor, 1e6))
        elif format_ == "cpm" and read_length > 0:
            # for cpm the size_factor is total reads
            self.__plus__ = np.divide(self.__plus__, np.divide(size_factor, 1e6))
            if self.__minus__ is not None:
                self.__minus__ = np.divide(self.__minus__, np.divide(size_factor, 1e6))
        elif size_factor is not None and size_factor > 0 and format_ == "atac":
            self.__plus__ = np.divide(self.__plus__, size_factor)  # * 100
            if self.__minus__ is not None:
                self.__minus__ = np.divide(self.__minus__, size_factor)


if __name__ == '__main__':
    pass

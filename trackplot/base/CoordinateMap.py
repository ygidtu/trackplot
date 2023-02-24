#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/11 2:55 PM

u"""
Convert aa position into genomic coordinate, Ran zhou.

"""
from itertools import islice
from typing import Optional

import numpy as np


class Coordinate(object):
    u"""
    A Coordinate object for genomic regions.
    """

    __slots__ = ["strand", "se"]

    def __init__(self, coordinate_list: list, strand: str = "*"):
        u"""
        Set genomic coordinates
        :param coordinate_list: a nested tuple of list
        :param strand: a strands of given coordinate object
        """
        self.strand = strand
        self.se = self.__fmt_exons__(coordinate_list)
        assert len(self.location_list) == len(set(self.location_list)), \
            f"Overlapped regions were found in {self.se}"

    @staticmethod
    def __fmt_exons__(coordinate_list: list) -> list:
        u"""
        Format and sort exon list
        :param coordinate_list: a nested tuple of list. like [('5','6'),('1','4')]
        :return: a nested tuple of list. like [(1,4), (5,6)]
        """
        # sorting coordinate based on first of location.
        formatted_coordinate_list = list(map(lambda x: tuple(map(int, x)), coordinate_list))
        formatted_coordinate_list.sort(key=lambda x: x[0])
        return formatted_coordinate_list

    @staticmethod
    def __get_s_or_e__(coordinate_list: list, index: int) -> list:
        u"""
        Get start or end site for each given coordinates
        :param coordinate_list: a nested tuple of list, like [(1,3),(5,6)]
        :param index: the index of tuple, 0 for the left site and 1 for the right end site.
        :return: a list which contained the start or end sites.
        """
        return list(map(lambda x: x[index], coordinate_list))

    @property
    def start(self):
        return self.__get_s_or_e__(coordinate_list=self.se, index=0)

    @property
    def end(self):
        return self.__get_s_or_e__(coordinate_list=self.se, index=1)

    @property
    def introns(self):
        u"""
        Set intronic regions for each coordinate object
        :return: a nested tuple of list which contained intronic coordinates.
        """
        if len(self.se) == 1:
            return None
        else:
            introns_list = []
            for left_exon, right_exon in self.__slide_window__(self.se, num_of_chunk=2):
                introns_list.append(tuple(
                    [left_exon[1] + 1,
                     right_exon[0] - 1]
                ))
            return introns_list

    @staticmethod
    def __slide_window__(nested_list: list, num_of_chunk: int):
        u"""
        A sliding window to slice the given list
        :param nested_list: a nested tuple of list, like [(1,3),(5,6)]
        :param num_of_chunk:  num of element for each chunk.
        :return: 
        """
        nested_list = iter(nested_list)
        chunked_list = list(islice(nested_list, num_of_chunk))
        if len(chunked_list) == num_of_chunk:
            yield chunked_list
        for elem in nested_list:
            result = chunked_list[1:] + list((elem,))
            yield result

    @classmethod
    def __flatten__(cls, nested_list: list):
        u"""
        Flatten the nested list
        :param nested_list: a nested tuple of list, like [(1,3),(5,6)]
        :return:
        """
        for sub_list in nested_list:
            if hasattr(sub_list, "__iter__") and not isinstance(sub_list, str):
                for sub_el in cls.__flatten__(sub_list):
                    yield sub_el
            else:
                yield sub_list

    @property
    def location_list(self) -> np.ndarray:
        u"""
        Get a list which contained all position.
        For example, an exon list `[(10, 15), (20,25)]` as input,
        [10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25] will return
        :return: a list which contained all position
        """
        # Add 1 offset because of 0-based coordinate
        position_list = list(self.__flatten__(list(map(
            lambda x: range(x[0], x[1] + 1), self.se
        ))))

        if self.strand == '-':
            return np.array(position_list[::-1])

        return np.array(position_list)

    @property
    def pep_index(self) -> np.ndarray:
        u"""
        Relative position of pep coordinates
        :return: 
        """
        return np.array(list(map(lambda x: int(x / 3), range(len(self.location_list)))))

    @property
    def cds_index(self) -> np.ndarray:
        u"""
        Relative position of cds coordinates
        :return:
        """
        return np.array(list(range(len(self.location_list))))

    @staticmethod
    def __group_consecutive_value__(location_list: np.ndarray, strand: str) -> list:
        u"""
        group the consecutive value into a list

        :param location_list: a list of location site
        :param strand: the strand of the current list to group
        :return:
        """
        offset = -1 if strand == "-" else 1
        group_ids = np.concatenate(([0], (np.diff(location_list) != offset).cumsum()))

        grouped_truncated_location_array = \
            np.split(
                location_list,
                np.unique(group_ids, return_counts=True)[1].cumsum().tolist()
            )

        return grouped_truncated_location_array

    @classmethod
    def init_from_location_list(cls, truncated_location_array: np.ndarray, strand: str):
        u"""
        init class based on the list of location
        :param truncated_location_array: truncated location array
        :param strand: the strand of the current list
        :return:
        """
        __coordinate_list = []
        for sub_array in cls.__group_consecutive_value__(truncated_location_array, strand):
            if len(sub_array) == 0:
                continue

            __coordinate_list.append(
                tuple(
                    [
                        min(sub_array),
                        max(sub_array)
                    ]
                )
            )
        return cls(__coordinate_list, strand)


class CoordinateMapper(Coordinate):
    u"""
    Convert positions between CDS and protein coordinates.
    TODO: Add cds to pep coordinate?
    """

    def __init__(self, coordinates_list, strand: str):
        u"""
        Set genomic coordinates to be used for mapping
        :param coordinates_list: a nested tuple of list
        :param strand: a strands of given coordinate object
        """
        super().__init__(coordinate_list=coordinates_list, strand=strand)

    def pep_to_cds(self, pep_start: int, pep_end: Optional[int] = None):
        u"""
        Convert pep position into genomic position
        :param pep_start: the start position of pep
        :param pep_end: the end position of pep, if None, the start site is equal to end site
        :return:
        """
        # change 1-based pep coordinate into 0-based coordinate
        pep_start = pep_start - 1
        pep_end = pep_start if pep_end is None else pep_end - 1

        start_ind = np.where(self.pep_index == pep_start)[0]
        end_ind = np.where(self.pep_index == pep_end)[0]

        cds_left_index, cds_right_index = start_ind[0], end_ind[-1] + 1

        return self.init_from_location_list(
            self.location_list[cds_left_index:cds_right_index],
            self.strand)


if __name__ == '__main__':
    pass

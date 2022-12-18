#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""


class Junction(object):
    u"""
    Created by ygidtu at 2018.12.19

    This is used to collect information of single junction
    And provide relative position comparison
    """

    __slots__ = ["chromosome", "start", "end", "strand"]

    def __init__(self, chromosome, start, end, strand: str = "+"):
        u"""
        init this class
        :param chromosome: the chromosome name of the given junction
        :param start: the start site of the given junction
        :param end: the end site of the given junction
        :param strand: the strand of the given junction.
        """
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

        if self.end <= self.start:
            raise ValueError(f"End site({start}) should bigger than start site({end})")

    @property
    def length(self):
        u"""
        :return: int, the length of this junction
        """
        return self.end - self.start

    @classmethod
    def create_junction(cls, string):
        u"""
        create Junction from chr1:1-100:+
        :param string: str, chr1:1-100:+ format or chr1:1-100 also work
        :return:
        """
        string = string.split(":")

        chromosome = string[0]
        start, end = string[1].split("-")

        strand = "+"
        if len(string) > 2:
            strand = string[-1]

        return cls(chromosome=chromosome, start=start, end=end, strand=strand)

    def __hash__(self):
        u"""
        generate hash
        :return:
        """
        return hash((self.chromosome, self.start, self.end, self.strand))

    def __str__(self):
        u"""
        convert junctions to string
        :return:
        """
        return "{chrom}:{start}-{end}".format(
            **{
                "chrom": self.chromosome,
                "start": self.start,
                "end": self.end
            }
        )

    def __gt__(self, other):
        u"""
        greater than
        compare two junction by length
        :param other:
        :return:
        """
        return self.length > other.length

    def __lt__(self, other):
        u"""
        less than
        compare two junction by length
        :param other:a
        :return:
        """
        return self.length < other.length

    def __eq__(self, other):
        u"""
        same length
        :param other:
        :return:
        """
        return self.length == other.length

    def is_overlap(self, other):
        u"""
        whether any overlap with another Junction or GenomicLoci
        :param other:
        :return:
        """

        if self.chromosome != other.chromosome:
            return False

        return self.start < other.end and self.end > other.start

    def is_upstream(self, other):
        u"""
        whether this junction is upstream of other
        :param other:
        :return:
        """
        assert isinstance(other, Junction), "Input should be Junction class"

        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome

        return self.end < self.start

    def is_downstream(self, other):
        u"""
        whether this junction is downstream of other
        :param other:
        :return:
        """
        assert isinstance(other, Junction), "Input should be Junction class"

        if self.chromosome != other.chromosome:
            return self.chromosome > other.chromosome

        return self.start > other.end


if __name__ == '__main__':
    pass

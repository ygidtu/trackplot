#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2022.05.11

The scripts contain the object to handle all the pysam, pybigwig related reading
"""
import os

import numpy as np
import pyBigWig
import pysam

from conf.logger import logger
from sashimi.base.GenomicLoci import GenomicLoci


def __opposite_strand__(strand: str) -> str:
    u"""
    replace to opposite of current strand information.
    :param strand: strand, one of ['+', '-', '*']
    :return: an opposite of strand
    """
    assert strand in ["+", "-", "*"], "Unknown strand information was found."
    if strand == "*":
        return "*"
    elif strand == "+":
        return "-"
    else:
        return "+"


def __get_strand__(read: pysam.AlignedSegment, library: str) -> str:
    u"""
    Determine the strand of for each read.
    :param read: a pysam.AlignedSegment from the bam file.
    :param library: the method for preparing of the library,
    value should be one of ["fr-firststrand", "fr-secondstrand", "fr-unstrand"]
    :return:
    """
    assert library in ["fr-firststrand", "fr-secondstrand", "fr-unstrand"], "Can't recognize the definition of library."

    # return '+' strand for all unstrand library.
    if library == "fr-unstrand":
        return "+"

    # Only guessing the strand based on fr-secondstrand rule.
    if read.is_paired:
        if read.is_read1 and read.is_reverse:
            current_strand = "+"
        elif read.is_read2 and not read.is_reverse:
            current_strand = "+"
        else:
            current_strand = "-"
    else:
        current_strand = "+" if read.is_reverse else "-"

    if library == "fr-secondstrand":
        return current_strand
    else:
        return __opposite_strand__(current_strand)


class Reader(object):

    @classmethod
    def read_bam(cls, path: str, region: GenomicLoci, library: str = "fr-unstrand"):
        chrom = region.chromosome
        if not os.path.exists(path + ".bai"):
            logger.info(f"try to create index for {path}")
            pysam.index(path)

        with pysam.AlignmentFile(path, 'rb') as bam_file:
            try:
                relevant_reads = bam_file.fetch(reference=chrom, start=region.start, end=region.end)
            except ValueError as err:
                logger.warning(err)

                if chrom.startswith("chr"):
                    logger.info("try without chr")
                    chrom = chrom.replace("chr", "")
                else:
                    logger.info("try with chr")
                    chrom = "chr{}".format(chrom)
                relevant_reads = bam_file.fetch(reference=chrom, start=region.start, end=region.end)

            for read in relevant_reads:
                yield read, __get_strand__(read, library=library)

    @classmethod
    def read_gtf(cls, path: str, region: GenomicLoci, bed: bool = False):
        with pysam.TabixFile(path) as r:
            try:
                iter_ = r.fetch(
                    region.chromosome,
                    region.start,
                    region.end,
                    parser=pysam.asGTF() if not bed else pysam.asBed())

            except ValueError:
                try:
                    if not region.chromosome.startswith("chr"):
                        logger.info("Guess need 'chr'")
                        iter_ = r.fetch(
                            "chr" + region.chromosome,
                            region.start, region.end,
                            parser=pysam.asGTF() if not bed else pysam.asBed()
                        )
                    else:
                        logger.info("Guess 'chr' is redundant")
                        iter_ = r.fetch(
                            region.chromosome.replace("chr", ""),
                            region.start,
                            region.end,
                            parser=pysam.asGTF() if not bed else pysam.asBed()
                        )
                except ValueError as err:
                    logger.warn(
                        "please check the input region and {} files".format("gtf" if not bed else "bed")
                    )
                    logger.error(err)
                    raise err

            for record in iter_:
                yield record

    @classmethod
    def read_bigwig(cls, path: str, region: GenomicLoci) -> np.array:
        with pyBigWig.open(path) as r:
            try:
                return r.values(region.chromosome, region.start, region.end + 1)
            except RuntimeError as e:
                logger.warning(e)

                logger.info("may be caused by the mismatch of chromosome")
                if region.chromosome.startswith("chr"):
                    return r.values(region.chromosome.replace("chr", ""), region.start, region.end + 1)
                else:
                    return r.values("chr" + region.chromosome, region.start, region.end + 1)


if __name__ == '__main__':
    pass

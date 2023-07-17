#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2022.05.11

The scripts contain the object to handle all the pysam, pybigwig related reading
"""
import os
import sys
from typing import Optional

import numpy as np
import pysam
from loguru import logger

from trackplot.base.GenomicLoci import GenomicLoci


def __opposite_strand__(strand: str) -> str:
    u"""
    replace to opposite of current strand information
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
        value should be one of [frf: "fr-firststrand", frs: "fr-secondstrand", fru: "fr-unstrand"]
    :return:
    """
    assert library in ["frf", "frs", "fru"], "Can't recognize the definition of library."

    # return '+' strand for all unstrand library.
    if library == "fru":
        return "+"

    # Only guessing the strand based on fr-firststrand rule.
    if read.is_paired:
        if read.is_read1 and read.is_reverse:
            current_strand = "+"
        elif read.is_read2 and not read.is_reverse:
            current_strand = "+"
        else:
            current_strand = "-"
    else:
        current_strand = "+" if read.is_reverse else "-"
        """
        TODO, add a more elegant way to solve some single-end strand specific library.
        """
        if library == "frf":
            # here for 3' tag rna seq
            current_strand = __opposite_strand__(current_strand)

    if library == "frf":
        return current_strand
    else:
        return __opposite_strand__(current_strand)


class Reader(object):

    @classmethod
    def __modify_chrom__(cls, region: GenomicLoci, reader, parser=None):
        if not region.chromosome.startswith("chr"):
            logger.info("Guess need 'chr'")

            if parser:
                iter_ = reader.fetch(
                    "chr" + region.chromosome if region.chromosome != "MT" else "chrM",
                    region.start, region.end, parser=parser
                )
            else:
                iter_ = reader.fetch(
                    "chr" + region.chromosome if region.chromosome != "MT" else "chrM",
                    region.start, region.end
                )

        else:
            logger.info("Guess 'chr' is redundant")
            if parser:
                iter_ = reader.fetch(
                    region.chromosome.replace("chr", "") if region.chromosome != "chrM" else "MT",
                    region.start, region.end, parser=parser
                )
            else:
                iter_ = reader.fetch(
                    region.chromosome.replace("chr", "") if region.chromosome != "chrM" else "MT",
                    region.start, region.end
                )

        return iter_

    @classmethod
    def read_bam(cls, path: str, region: GenomicLoci, library: str = "fru"):
        chrom = region.chromosome
        if not os.path.exists(path + ".bai"):
            logger.info(f"try to create index for {path}")
            pysam.index(path)

        with pysam.AlignmentFile(path, 'rb') as bam_file:
            try:
                relevant_reads = bam_file.fetch(reference=chrom, start=region.start, end=region.end)
            except ValueError as err:
                # logger.debug(err)
                relevant_reads = cls.__modify_chrom__(region, bam_file)

            for read in relevant_reads:
                if read.is_qcfail or read.is_unmapped or read.is_duplicate:
                    continue

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
                    iter_ = cls.__modify_chrom__(region, r, parser=pysam.asGTF() if not bed else pysam.asBed())
                except ValueError as err:
                    logger.debug("please check the input region and gtf files")
                    logger.error(err)
                    raise err

            for record in iter_:
                yield record

    @classmethod
    def read_bigwig(cls, path: str, region: GenomicLoci) -> np.array:
        try:
            import pyBigWig
        except ImportError as err:
            logger.error(f"Please install pyBigWig properly to enable bigWig support: {err}")
            sys.exit(1)

        with pyBigWig.open(path) as r:
            try:
                return r.values(region.chromosome, region.start, region.end + 1)
            except RuntimeError as e:
                logger.debug(e)

                logger.info("may be caused by the mismatch of chromosome")
                if region.chromosome.startswith("chr"):
                    return r.values(region.chromosome.replace("chr", ""), region.start, region.end + 1)
                else:
                    return r.values("chr" + region.chromosome, region.start, region.end + 1)

    @classmethod
    def read_bigbed(cls, path: str, region: GenomicLoci):
        try:
            import pyBigWig
        except ImportError as err:
            logger.error(f"Please install pyBigWig properly to enable bigWig support: {err}")
            sys.exit(1)

        with pyBigWig.open(path) as r:
            if not r.isBigBed():
                logger.debug(f"{path} don't look like bigbed file.")
                yield from []
            else:
                try:
                    iter_ = r.entries(region.chromosome, region.start, region.end + 1)
                except RuntimeError as e:
                    logger.debug(e)

                    logger.info("may be caused by the mismatch of chromosome")
                    if region.chromosome.startswith("chr"):
                        iter_ = r.entries(
                            region.chromosome.replace("chr", "") if region.chromosome != "chrM" else "MT",
                            region.start,
                            region.end + 1
                        )
                    else:
                        iter_ = r.entries(
                            "chr" + region.chromosome if region.chromosome != "MT" else "chrM",
                            region.start,
                            region.end + 1
                        )
                if not iter_:
                    iter_ = []
                yield from iter_

    @classmethod
    def read_depth(cls, path: str, region: Optional[GenomicLoci] = None):
        if not os.path.exists(path + ".tbi"):
            logger.info(f"create tbi index for {path}")
            pysam.tabix_index(
                path, seq_col=0,
                start_col=1, end_col=1,
                force=True, keep_original=True
            )

        with pysam.TabixFile(path) as r:
            if not region:
                iter_ = r.fetch()
            else:
                try:
                    iter_ = r.fetch(region.chromosome, region.start, region.end, parser=pysam.asTuple())
                except ValueError:
                    try:
                        iter_ = cls.__modify_chrom__(region, r, parser=pysam.asTuple())
                    except ValueError as err:
                        logger.debug("please check the input region and gtf files")
                        logger.error(err)
                        raise err

            for record in iter_:
                yield record

    @classmethod
    def read_hic(cls, path: str, region: GenomicLoci):
        try:
            from hicmatrix import HiCMatrix as hm
        except ImportError as err:
            logger.error(f"Please install pyBigWig properly to enable HiC support: {err}")
            sys.exit(1)

        os.environ['NUMEXPR_MAX_THREADS'] = '16'
        os.environ['NUMEXPR_NUM_THREADS'] = '8'
        hic = hm.hiCMatrix(path, f"{region.chromosome}:{region.start}-{region.end}")
        return hic

    @classmethod
    def total_reads_of_bam(cls, path: str):
        if not os.path.exists(path + ".bai"):
            logger.info(f"try to create index for {path}")
            pysam.index(path)

        total = 0
        with pysam.AlignmentFile(path, 'rb') as bam_file:
            for _ in bam_file:
                total += 1
        return total


if __name__ == '__main__':
    pass

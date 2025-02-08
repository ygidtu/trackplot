#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
This file contains the object to handle bam file related issues.

changelog:
    1. add library parameter for determining of read strand at 2022.4.28.

"""
import gzip
import os
from typing import Optional, Set

import numpy as np
import pysam
from loguru import logger

from trackplot.base.GenomicLoci import GenomicLoci
from trackplot.base.Junction import Junction
from trackplot.base.ReadDepth import ReadDepth
from trackplot.base.Readder import Reader
from trackplot.conf.config import NORMALIZATION
from trackplot.file.File import SingleCell  


class Bam(SingleCell):

    __slots__ = "title", "label", "library", "density_by_strand"

    def __init__(self,
                 path: str, label: str = "",
                 title: str = "", barcodes: Optional[Set[str]] = None,
                 barcode_tag: str = "CB", umi_tag: str = "UB",
                 library: str = "fru", density_by_strand: bool = False,
                 size_factor: Optional[int] = None):
        u"""
        init this object
        :param label: the left axis label
        :param title: the default title to show in the upper-right of density plot
        :param barcodes: the path to barcodes,
                default: ../filtered_feature_bc_matrix/barcodes.tsv.gz of bam file according to 10X Genomics
        :param barcode_tag: the cell barcode tag, default is CB according to 10X Genomics
        :param umi_tag: the UMI barcode tag, default is UB according to 10X Genomics
        :param library: library for determining of read strand.
        :param density_by_strand: whether to draw density plot in strand-specific manner.
        """
        super().__init__(path, barcodes, barcode_tag, umi_tag)
        self.title = title
        self.label = label if label else os.path.basename(path).replace(".bam", "")
        self.library = library
        self.density_by_strand = density_by_strand
        self.size_factor = size_factor

    @classmethod
    def create(cls,
               path: str,
               label: str = "",
               title: str = "",
               barcodes: Optional[Set[str]] = None,
               barcode_tag: str = "CB",
               umi_tag: str = "UB",
               library: str = "fru",
               density_by_strand: bool = False,
               size_factor: Optional[int] = None
               ):
        u"""

        :param path: the path to bam file
        :param label: the left axis label
        :param title: the default title to show in the upper-right of density plot
        :param barcodes: the path to barcodes,
                default: ../filtered_feature_bc_matrix/barcodes.tsv.gz of bam file according to 10X Genomics
        :param barcode_tag: the cell barcode tag, default is CB according to 10X Genomics
        :param umi_tag: the UMI barcode tag, default is UB according to 10X Genomics
        :param library: library for determining of read strand.
        :param density_by_strand:
        :param size_factor:
        :return:
        """

        if not os.path.exists(path + ".bai"):
            pysam.index(path)

        barcode = barcodes
        path = os.path.abspath(path)
        if not barcodes:
            barcode = set()
            barcodes = os.path.join(os.path.dirname(path), "filtered_feature_bc_matrix/barcodes.tsv.gz")

            if os.path.exists(barcodes):
                with gzip.open(barcodes, "rt") as r:
                    for line in r:
                        barcode.add(line.strip())
        return cls(
            path=path,
            label=label,
            title=title,
            barcodes=barcode,
            barcode_tag=barcode_tag,
            umi_tag=umi_tag,
            library=library,
            density_by_strand=density_by_strand,
            size_factor=size_factor
        )

    def __hash__(self):
        return hash(self.label)

    def __str__(self) -> str:

        temp = []

        for x in [self.title, self.label, self.path]:
            if x is None or x == "":
                x = "None"
            temp.append(str(x))

        return "\t".join(temp)

    def to_csv(self) -> str:
        temp = []

        for x in [self.title, self.label, self.path]:
            if x is None or x == "":
                x = "None"
            if isinstance(x, list):
                x = ";".join(x)
            temp.append(str(x))

        return ",".join(temp)

    def load(self,
             region: GenomicLoci,
             threshold: int = 0,
             reads1: Optional[bool] = None,
             required_strand: Optional[str] = None,
             normalize_format: Optional[str] = None,
             **kwargs
             ):
        """
            determine_depth determines the coverage at each base between start_coord and end_coord, inclusive.

            bam_file_path is the path to the bam file used to \
            determine the depth and junctions on chromosome between start_coord and end_coord

        return values:
            depth_vector,
            which is a Numpy array which contains the coverage at each base position between start_coord and end_coord

            spanned_junctions, which is a dictionary containing the junctions supported by reads.
            The keys in spanned_junctions are the
                names of the junctions, with the format chromosome:lowerBasePosition-higherBasePosition
        :param region: GenomicLoci object including the region for calculating coverage
        :param threshold: minimums counts of the given splice junction for visualization
        :param reads1: None -> all reads, True -> only R1 kept; False -> only R2 kept
        :param required_strand: None -> all reads, else reads on specific strand
        :param normalize_format: None -> raw counts; others fpkm and cpm
        """
        self.region = region

        spanned_junctions = kwargs.get("junctions", {})
        included_junctions = kwargs.get("included_junctions", {})
        remove_duplicate_umi = kwargs.get("remove_duplicate_umi", False)
        spanned_junctions_plus, spanned_junctions_minus = {}, {}
        plus, minus = np.zeros(len(region), dtype=np.int32), np.zeros(len(region), dtype=np.int32)
        site_plus, site_minus = np.zeros(len(region), dtype=np.int32), np.zeros(len(region), dtype=np.int32)

        umis = {}
        read_lens = []
        scatac_alert = True
        try:
            for read, strand in Reader.read_bam(path=self.path, region=region, library=self.library):
                # make sure that the read can be used
                cigar_string = read.cigartuples

                # each read must have a cigar string
                if cigar_string is None:
                    continue

                # select R1 or R2
                if reads1 is True and not read.is_read1:
                    continue

                if reads1 is False and not read.is_read2:
                    continue

                if normalize_format != NORMALIZATION[0] and normalize_format is not None:
                    read_lens.append(read.query_alignment_length)

                # filter reads by 10x barcodes
                # @20220924, add `not` before has_barcode and skip these reads without umi tag.
                if self.barcodes:
                    if not read.has_tag(self.barcode_tag) or not self.has_barcode(read.get_tag(self.barcode_tag)):
                        continue

                    if remove_duplicate_umi:
                        barcode = read.get_tag(self.barcode_tag)
                        if barcode not in umis.keys():
                            umis[barcode] = {}

                        # filter reads with duplicate umi by barcode
                        if read.has_tag(self.umi_tag):
                            umi = read.get_tag(self.umi_tag)
                        else:
                            if scatac_alert:
                                logger.debug(f"Is {self.path} an scATAC-seq bam? There is no {self.umi_tag}")
                                scatac_alert = False
                            umi = read.query_name

                        if umi in umis[barcode].keys() and umis[barcode][umi] != hash(read.query_name):
                            continue

                        if len(umis[barcode]) == 0:
                            umis[barcode][umi] = hash(read.query_name)
                        # else:
                        #     # There is no umi tag in atacdata, and we add it
                        #     continue

                start = read.reference_start
                if required_strand and strand != required_strand:
                    continue

                """
                M	BAM_CMATCH	0
                I	BAM_CINS	1
                D	BAM_CDEL	2
                N	BAM_CREF_SKIP	3
                S	BAM_CSOFT_CLIP	4
                H	BAM_CHARD_CLIP	5
                P	BAM_CPAD	6
                =	BAM_CEQUAL	7
                X	BAM_CDIFF	8
                B	BAM_CBACK	9
                """
                for cigar, length in cigar_string:
                    cur_start = start + 1
                    cur_end = start + length + 1

                    if cigar == 0:  # M
                        for i in range(length):
                            if region.start <= start + i + 1 <= region.end:
                                try:
                                    if strand == "+":
                                        plus[start + i + 1 - region.start] += 1
                                    elif strand == "-":
                                        minus[start + i + 1 - region.start] += 1
                                    else:
                                        pass
                                except IndexError as err:
                                    logger.info(region)
                                    logger.info(cigar_string)
                                    logger.info(start, i)
                                    exit(err)

                    # remove the deletion.
                    if cigar not in (1, 4, 5):  # I, S, H
                        start += length

                    if cigar == 3 and not kwargs.get("only_customized_junction"):  # N
                        try:
                            junction_name = Junction(region.chromosome, cur_start, cur_end, strand)

                            if junction_name not in spanned_junctions:
                                spanned_junctions[junction_name] = 0

                            spanned_junctions[junction_name] = spanned_junctions[junction_name] + 1
                        except ValueError as err:
                            logger.debug(err)
                            continue
                start = read.reference_start + 1 if read.reference_start + 1 > region.start else region.start
                end = read.reference_end + 1 if read.reference_end + 1 < region.end else region.end
                if strand == "+" and 0 <= start - region.start < len(plus):
                    site_plus[end - region.start] += 1
                elif strand == "-" and 0 <= end - region.start < len(minus):
                    site_minus[start - region.start] += 1

            for k, v in spanned_junctions.items():
                kept = v >= threshold
                
                # if the number of junctiosn is lower than threshold, then skip
                if not kept:
                    continue
                
                # if included_junctions is provided, then skip all junctions by default
                if included_junctions:
                    kept = False
                
                # check whether junctions should be kept
                if k.str(with_strand = True) in included_junctions:
                    logger.debug(f"{str(k)} is included")
                    kept = True
                elif k.str(with_strand = False) in included_junctions:
                    logger.debug(f"{str(k)} is included, but strand is ignored")
                    kept = True

                if not kept:
                    logger.debug(f"{str(k)} is not included")
                else:
                    if k.strand == "+":
                        spanned_junctions_plus[k] = 1 + spanned_junctions_plus.get(k, v)
                    elif k.strand == "-":
                        spanned_junctions_minus[k] = -1 + spanned_junctions_minus.get(k, v)
                        
        except IOError as err:
            logger.error('There is no .bam file at {0}'.format(self.path))
            logger.error(err)
        except ValueError as err:
            logger.error(self.path)
            logger.error(err)

        self.data = ReadDepth(
            plus if self.density_by_strand else plus + minus,
            site_plus=site_plus,
            site_minus=site_minus,
            minus=minus if self.density_by_strand else None,
            junction_dict_plus=spanned_junctions_plus,
            junction_dict_minus=spanned_junctions_minus,
            strand_aware=False if self.library == "fru" else True)

        if normalize_format != NORMALIZATION[0] and normalize_format is not None:
            if self.size_factor is None:
                logger.info(f"Counting total number of reads: {self.path}")
                self.size_factor = Reader.total_reads_of_bam(self.path)
            else:
                logger.info(f"Using given total number of reads [{self.size_factor}] for {self.path}")
            self.data.normalize(size_factor=self.size_factor, format_=normalize_format, read_length=np.mean(read_lens))
        return self


if __name__ == '__main__':
    pass

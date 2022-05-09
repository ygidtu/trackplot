#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06

changelog:
    1. add library parameter for determining of read strand at 2022.4.28.

"""

from typing import Optional, List

import numpy as np
import pysam
from scipy.stats import zscore

from conf.logger import logger
from sashimi.file.BamInfo import BamInfo
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.Junction import Junction
from sashimi.base.Transcript import Transcript


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


class ReadDepth(GenomicLoci):
    u"""
    Migrated from SplicePlot ReadDepth class

    add a parent class to handle all the position comparison
    """

    def __init__(self,
                 chromosome,
                 start,
                 end,
                 wiggle,
                 junctions_dict,
                 reads=None,
                 plus=None,
                 minus=None,
                 library="fr-unstrand"):
        u"""
        init this class

        :param chromosome: the chromosome id of interesting region.
        :param start: the left site of genomic coordinate.
        :param end: the right site of genomic coordinate.
        :param wiggle: a numpy.ndarray object represented the whole read coverage.
        :param junctions_dict: a dict represented the coordinate of each intron as well as frequency.
        :param reads:
        :param plus:
        :param minus:
        :param library: the method for preparing of the library.
        """
        super().__init__(chromosome, start, end, "+")

        if wiggle is not None:
            assert chromosome is None or self.length + 1 == len(wiggle), "Wiggle length don't correspond to input range"

        self.wiggle = wiggle
        self.junctions_dict = junctions_dict
        self.max = max(self.wiggle)
        self.sequence = None
        self.chromosome = chromosome
        self.start = start
        self.__reads__ = reads if reads is not None else []
        self.plus = plus
        self.minus = minus * -1 if minus is not None else minus
        self.library = library

    @classmethod
    def determine_depth(
        cls,
        bam: BamInfo,
        chrom: str,
        start_coord: int,
        end_coord: int,
        threshold: int,
        threshold_of_reads: int,
        log,
        reads1: Optional[bool] = None,
        barcode_tag: str = "CB",
        umi_tag: Optional[str] = None,
        required_strand: Optional[str] = None,
        stack: bool = False,
        library: str = "fr-unstrand",
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

        :param chrom:
        :param barcode_tag:
        :param threshold_of_reads:
        :param threshold:
        :param end_coord:
        :param start_coord:
        :param bam:
        :param log: whether to use log transformed number of reads in sashimi
        :param reads1: None -> all reads, True -> only R1 kept; False -> only R2 kept
        :param required_strand: None -> all reads, else reads on specific strand
        :param stack: whether to kept reads info for stack plot
        :param library: the method for preparing of the library.
        """
        reads = {}
        filtered_reads = set()
        filtered_junctions = {}
        depth_vector = np.zeros(end_coord - start_coord + 1, dtype='f')
        spanned_junctions = {}
        plus, minus = np.zeros(end_coord - start_coord + 1, dtype="f"), np.zeros(end_coord - start_coord + 1, dtype="f")

        for bam_file_path in bam.path:
            try:
                with pysam.AlignmentFile(bam_file_path, 'rb') as bam_file:
                    try:
                        relevant_reads = bam_file.fetch(reference=chrom, start=start_coord, end=end_coord)
                    except ValueError as err:
                        logger.warning(err)
                        err = str(err)

                        if "without index" in err:
                            logger.info(f"try to create index for {bam_file_path}")
                            pysam.index(bam_file_path)
                            relevant_reads = bam_file.fetch(reference=chrom, start=start_coord, end=end_coord)
                        else:
                            if chrom.startswith("chr"):
                                logger.info("try without chr")
                                chrom = chrom.replace("chr", "")
                            else:
                                logger.info("try with chr")
                                chrom = "chr{}".format(chrom)
                            relevant_reads = bam_file.fetch(reference=chrom, start=start_coord, end=end_coord)

                    # tqdm()
                    for read in relevant_reads:
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

                        # filter reads by 10x barcodes
                        if not bam.empty_barcode():
                            if not read.has_tag(barcode_tag) or not bam.has_barcode(read.get_tag(barcode_tag)):
                                continue

                        start = read.reference_start
                        strand = __get_strand__(read=read, library=library)

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
                        exons_in_read = []
                        for cigar, length in cigar_string:
                            cur_start = start + 1
                            cur_end = start + length + 1

                            if cigar == 0:  # M
                                for i in range(length):
                                    if start_coord <= start + i + 1 <= end_coord:
                                        try:
                                            depth_vector[start + i + 1 - start_coord] += 1
                                        except IndexError as err:
                                            logger.info(start_coord, end_coord)
                                            logger.info(cigar_string)
                                            logger.info(start, i)
                                            exit(err)

                                if cur_start < end_coord and cur_end > start_coord:
                                    exons_in_read.append(GenomicLoci(
                                        chromosome=read.reference_name,
                                        start=cur_start if cur_start > start_coord else start_coord,
                                        end=cur_end if cur_end <= end_coord else end_coord,
                                        strand="+",
                                    ))

                            if cigar not in (1, 2, 4, 5):  # I, D, S, H
                                start += length

                            if cigar == 3:  # N
                                try:
                                    junction_name = Junction(chrom, cur_start, cur_end)

                                    if junction_name not in spanned_junctions:
                                        spanned_junctions[junction_name] = 0

                                    spanned_junctions[junction_name] = spanned_junctions[junction_name] + 1
                                except ValueError as err:
                                    logger.warning(err)
                                    continue

                        t = Transcript(
                            chromosome=read.reference_name,
                            start=read.reference_start + 1 if read.reference_start + 1 > start_coord else start_coord,
                            end=read.reference_end + 1 if read.reference_end + 1 < end_coord else end_coord,
                            strand=strand,
                            exons=exons_in_read,
                            is_reads=True
                        )

                        if stack:
                            reads[t] = reads.get(t, 0) + 1

                        if strand == "+" and read.reference_start >= start_coord:
                            plus[t.start - start_coord] += 1
                        elif strand == "-" and read.reference_end <= end_coord:
                            minus[t.end - start_coord] += 1

                for k, v in spanned_junctions.items():
                    if v >= threshold:
                        filtered_junctions[k] = v

                if log == 10:
                    depth_vector = np.log10(depth_vector + 1)
                elif log == 2:
                    depth_vector = np.log2(depth_vector + 1)
                elif log == "zscore":
                    depth_vector = zscore(depth_vector)

                for k, v in reads.items():
                    if v >= threshold_of_reads:
                        filtered_reads.add(k)
            except IOError as err:
                logger.error('There is no .bam file at {0}'.format(bam_file_path))
                logger.error(err)
            except ValueError as err:
                logger.error(bam_file_path)
                logger.error(err)

        if bam.show_mean:
            depth_vector = depth_vector / len(bam.path)
            plus = plus / len(bam.path)
            minus = minus / len(bam.path)

        return cls(
            chromosome=chrom,
            start=start_coord,
            end=end_coord,
            wiggle=depth_vector,
            junctions_dict=filtered_junctions,
            reads=filtered_reads,
            plus=plus,
            minus=minus
        )

    @classmethod
    def __tabix_iter__(cls, tbx, chrom: str, start_coord: int, end_coord: int):
        try:
            r = tbx.fetch(chrom, start_coord, end_coord, parser=pysam.asTuple())
        except ValueError as err:
            logger.warning(err)

            if chrom.startswith("chr"):
                logger.info("try without chr")
                chrom = chrom.replace("chr", "")
            else:
                logger.info("try with chr")
                chrom = "chr{}".format(chrom)

            r = tbx.fetch(chrom, start_coord, end_coord, parser=pysam.asTuple())
        return r

    @classmethod
    def determine_depth_by_fragments(
        cls, bam: BamInfo,
        chrom: str,
        start_coord: int,
        end_coord: int,
        strand: str,
        strandless: bool,
        log,
    ):
        filtered_reads = set()
        filtered_junctions = {}
        depth_vector = np.zeros(end_coord - start_coord + 1, dtype='f')

        for bam_file_path in bam.path:
            with pysam.TabixFile(bam_file_path) as tbx:
                for _, start, end, barcode, count in cls.__tabix_iter__(tbx, chrom, start_coord, end_coord):
                    # filter reads by 10x barcodes
                    if not bam.empty_barcode():
                        if not bam.has_barcode(barcode):
                            continue

                    start, end, count = int(start),  int(end), int(count)
                    for i in range(max(start, start_coord), min(end, end_coord)):
                        depth_vector[i - start_coord] += count

        if bam.show_mean:
            depth_vector = depth_vector / len(bam.path)

        if log == 10:
            depth_vector = np.log10(depth_vector + 1)
        elif log == 2:
            depth_vector = np.log2(depth_vector + 1)
        elif log == "zscore":
            depth_vector = zscore(depth_vector)

        return cls(
            chromosome=chrom,
            start=start_coord,
            end=end_coord,
            wiggle=depth_vector,
            junctions_dict=filtered_junctions,
            reads=filtered_reads,
            plus=depth_vector if strandless or strand == "+" else None,
            minus=depth_vector if strandless or strand == "-" else None,
        )

    @classmethod
    def determine_depth_by_depths(
            cls, bam: BamInfo,
            chrom: str,
            start_coord: int,
            end_coord: int,
            strand: str,
            strandless: bool,
            log,
            groups: List[int] = []
    ):
        u"""
        depth file generated by samtools depth

        samtools depth *.bam | bgzip > depth.gz && tabix indexed tabix -s 1 -b 2 -e 3 depth.gz
        :param bam:
        :param chrom:
        :param start_coord:
        :param end_coord:
        :param strand:
        :param strandless:
        :param log:
        :param groups: Used to query specific data from depth file
        :return:
        """
        depth_vector = np.zeros(end_coord - start_coord + 1, dtype='f')
        for bam_file_path in bam.path:
            with pysam.TabixFile(bam_file_path) as tbx:
                for vals in cls.__tabix_iter__(tbx, chrom, start_coord, end_coord):
                    chrom, site = vals[0], int(vals[1])

                    if not groups:
                        count = sum([int(x) for x in vals[3:]])
                    else:
                        count = 0
                        for g in groups:
                            try:
                                count += int(vals[g])
                            except IndexError as e:
                                logger.warning(f"the {g} is out of bound [{len(vals)}]: {e}")

                    depth_vector[site - start_coord] += count

        if bam.show_mean:
            depth_vector /= len(bam.path)

        if log == 10:
            depth_vector = np.log10(depth_vector + 1)
        elif log == 2:
            depth_vector = np.log2(depth_vector + 1)
        elif log == "zscore":
            depth_vector = zscore(depth_vector)

        return cls(
            chromosome=chrom,
            start=start_coord,
            end=end_coord,
            wiggle=depth_vector,
            junctions_dict={},
            reads={},
            plus=depth_vector if strandless or strand == "+" else None,
            minus=depth_vector if strandless or strand == "-" else None,
        )

    @classmethod
    def create_depth(cls, data, splice_region, depth=10):
        u"""
        Create ReadDepth base on junction dict
        :param data: {junction in string: int}
        :param splice_region:
        :param depth:
        :return:
        """
        junctions_dict = {}
        for key, value in data.items():
            junctions_dict[Junction.create_junction(key)] = value

        depth_vector = np.zeros(
            splice_region.end - splice_region.start + 1,
            dtype='f'
        )

        for i in splice_region.exons:
            for j in range(i.start - splice_region.start + 1, i.end - splice_region.start + 2):
                depth_vector[j] = depth

        return cls(
            chromosome=splice_region.chromosome,
            start=splice_region.start,
            end=splice_region.end,
            wiggle=depth_vector,
            junctions_dict=junctions_dict
        )

    @classmethod
    def create_blank(cls):

        """
            create_blank creates an instance of ReadDepth where all the attributes are None
        """
        return cls("1", 0, 1, None, None)

    def is_invalid(self):
        """
            is_invalid determines whether any of the attributes are None
        """
        return self.junctions_dict is None or (self.chromosome == "1" and self.start == 0 and self.end == 1)

    def shrink(self, new_low, new_high):

        """
            shrink changes the boundaries of the ReadDepth object

            new_low is the new lower genomic coordinate boundary for the ReadDepth object
            new_high is the new upper genomic coordinate boundary for the ReadDepth object

            This method also changes self.wiggle and self.
            junctions_dict so that they only contain data between new_low and new_high

            return value:
                Nothing. Method changes the ReadDepth object
        """
        if new_low < self.start or new_high > self.end:
            raise Exception(f'New boundaries are not valid, old: {self.start}-{self.end}, new: {new_low}-{new_high}')

        # filter through junctions_dict to remove junctions which are no longer in the region
        new_junctions_dict = {}
        for key, value in self.junctions_dict.items():
            ss_low, ss_high = key.start, key.end

            if ss_low >= new_low and ss_high <= new_high:
                new_junctions_dict[key] = value

        self.junctions_dict = new_junctions_dict

        # replace the wiggle
        bottom_index = new_low - self.start
        top_index = bottom_index + (new_high - new_low)

        self.wiggle = self.wiggle[bottom_index:top_index + 1]

        # change the lower and upper bound (last step)
        self.start = new_low
        self.end = new_high

    def __add__(self, other):

        """
            __add__ allows two ReadDepth objects to be added together using the + symbol

            Both self and other must have the same low and high attributes

            return value:
                A new ReadDepth object containing the sum of the two original ReadDepth objects
        """

        if self.is_invalid():
            return other
        if other.is_invalid():
            return self

        assert self.chromosome == other.chromosome, 'Cannot add depths from different chromosomes'
        assert self.start == other.start and self.end == other.end, 'Cannot add depths with different start and end'
        new_wiggle = self.wiggle + other.wiggle

        return ReadDepth(
            self.chromosome,
            self.start,
            self.end,
            new_wiggle,
            self.add_customized_junctions(other)
        )

    def __str__(self):
        return '{0}:{1}-{2},{3},{4}'.format(
            self.chromosome,
            self.start,
            self.end,
            self.wiggle,
            self.junctions_dict
        )

    def divide_by_constant(self, constant):
        """
        divide_by_constant divides self.wiggle and self.junctions_dict by a constant value

        constant is a number

        return value:
            A new ReadDepth object containing the divided values. Method leaves the original ReadDepth object unchanged
        """
        new_wiggle = self.wiggle / constant

        new_junctions_dict = {}
        for key, value in list(self.junctions_dict.items()):
            new_junctions_dict[key] = value * 1.0 / constant

        return ReadDepth(
            self.chromosome,
            self.start,
            self.end,
            new_wiggle,
            new_junctions_dict
        )

    def filter_junctions_dict_for_event(self, splice_event_name):
        """
            filter_junctions_dict_for_event removes all entries frm junctions_dict that cannot possibly be
                involved in the alternative splicing event splice_event_name

            splice_event_name is the name of the alternative splicing event,
            in the format -> chr1:17055-17915,chr1:17055-17606,chr1:17055-17233,
                where the numbers represent the genomic coordinates of the splice sites

            return values:
                A new ReadDepth object containing only the relevant junctions
        """
        junction_names_list = splice_event_name.split(',')
        new_junctions_dict = {}
        for junction_name in junction_names_list:
            if junction_name in self.junctions_dict:
                new_junctions_dict[junction_name] = self.junctions_dict[junction_name]

        return ReadDepth(
            self.chromosome,
            self.start,
            self.end,
            self.wiggle,
            new_junctions_dict
        )

    def get_read_depth(self, genomic):
        u"""
        Created by ygidtu at 2018.12.25
        Extract part of data from this class
        :param genomic: a genomic class
        :return:
        """
        if genomic.start == self.start and genomic.end == self.end:
            return self

        targets = [x for x in range(genomic.start - self.start, genomic.end - self.start + 1)]
        wiggle = self.wiggle.take(targets)

        junctions_dict = {}
        for k, v in self.junctions_dict.items():
            if k.is_overlap(genomic):
                junctions_dict[k] = v

        return ReadDepth(
            chromosome=self.chromosome,
            start=genomic.start,
            end=genomic.end,
            wiggle=wiggle,
            junctions_dict=junctions_dict
        )

    def add_customized_junctions(self, other):
        u"""
        Add customized junctions to plot
        :param other:
        :return:
        """
        new_junctions_dict = {}

        for key, value in self.junctions_dict.items():
            if key in other.junctions_dict:
                new_junctions_dict[key] = value + other.junctions_dict[key]
            else:
                new_junctions_dict[key] = value

        for key, value in list(other.junctions_dict.items()):
            if key not in self.junctions_dict:
                new_junctions_dict[key] = value

        self.junctions_dict = new_junctions_dict
        return new_junctions_dict

    def __iter__(self):
        for idx, val in enumerate(self.wiggle):
            for i in range(int(val)):
                yield "{},{},{}".format(self.chromosome, self.start + idx, val)

    @property
    def reads(self):
        u"""
        convert self.__transcripts__ to list format

        Note: wrapper gtf proxy to list and dict format
            1. there no need to change the code of sashimi plot
            2. shrink the transcript range

        :return: [[{transcript: id, gene: id, exon: []}, {}, {}], [{}]]
        """
        return sorted(
            self.__reads__,
            key=lambda x: (x.start, x.end, len(x.exons)),
            reverse=True
        )

    @property
    def exon_starts(self):
        u"""
        API for extract all exon starts
        :return:
        """
        starts = []
        for i in self.reads:
            for j in i.exons:
                starts.append(j.start)
        return sorted(starts)

    @property
    def exon_ends(self):
        u"""
        API for extract all exon ends
        :return:
        """
        ends = []
        for i in self.reads:
            for j in i.exons:
                ends.append(j.end)
        return sorted(ends)


if __name__ == '__main__':
    pass

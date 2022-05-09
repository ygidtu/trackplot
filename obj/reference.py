#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
This file contains the reference file definition and related functions
"""
import gzip
import os
import re

import filetype
import pysam

from conf.logger import logger
from src.GenomicLoci import GenomicLoci
from src.Transcript import Transcript


class Reference(object):

    u"""
    The reference file, support gtf and gff format
    """

    def __init__(self, path: str):
        u"""
        init func
        :param path: path to input file
        """
        self.path = path
        self.transcripts = []

    def __add__(self, other):
        assert isinstance(other, Reference), "only Reference and Reference could be added"

        self.transcripts += other.transcripts
        self.transcripts = sorted(set(self.transcripts))

    @classmethod
    def create(cls, path: str, region: GenomicLoci):
        u"""
        create reference file object
        :param path: path to input file
        :param region: target region
        :return: Reference obj
        """
        assert os.path.exists(path), f"{path} not exists"

        obj = cls(path=cls.index_gtf(path))
        obj.__load__(region)
        return obj

    @staticmethod
    def is_gtf(infile):
        u"""
        check if input file is gtf
        :param infile: path to input file
        :return:
        """
        if infile is None:
            return False

        is_gtf = 0
        try:
            if filetype.guess_mime(infile) == "application/gzip":
                is_gtf += 10
                r = gzip.open(infile, "rt")
            else:
                r = open(infile)

            for line in r:
                if line.startswith("#"):
                    continue

                lines = re.split(r"\s+", line)

                if len(lines) < 8:
                    break

                if re.search(r"([\w-]+ \"[\w.\s\-%,:]+\";? ?)+", " ".join(lines[8:])):
                    is_gtf += 1

                break

            r.close()
            return is_gtf
        except TypeError as err:
            logger.error("failed to open %s", infile)
            exit(err)

    @classmethod
    def index_gtf(cls, input_gtf, sort_gtf=True, retry=0):
        u"""
        Check the tabix index of input gtf file

        Extract only exon tags and keep it clean
        :param input_gtf: path to input gtf file
        :param sort_gtf: Boolean value, whether to sort gtf file first
        :param retry: only try to sort gtf once
        :return path to compressed and indexed bgzipped gtf file
        """
        gtf = cls.is_gtf(input_gtf)
        assert gtf % 10 == 1, f"{input_gtf} seems not be gtf format"

        index = False
        if gtf // 10 > 0:
            output_gtf = input_gtf
        else:
            output_gtf = input_gtf + ".gz"

        sorted_gtf = re.sub(r"\.gtf$", "", input_gtf) + ".sorted.gtf"

        if not os.path.exists(output_gtf) or not os.path.exists(output_gtf + ".tbi"):
            logger.info(f"the {output_gtf} or tbi index not exists")
            index = True
        elif os.path.getctime(output_gtf) < os.path.getctime(output_gtf) or \
                os.path.getctime(output_gtf) < os.path.getctime(output_gtf):

            logger.info("the tbi index is older than the gtf file")
            index = True

        if not index:
            return output_gtf

        logger.info("Create index for %s", input_gtf)
        try:
            pysam.tabix_index(
                input_gtf,
                preset="gff",
                force=True,
                keep_original=True
            )
        except OSError as err:

            if re.search("could not open", str(err)):
                raise err

            logger.error(err)
            logger.error("Guess gtf needs to be sorted")

            data = []
            logger.info("Sorting %s" % input_gtf)

            old_input_gtf = input_gtf
            input_gtf = re.sub(r"\.gtf$", "", input_gtf) + ".sorted.gtf"
            output_gtf = input_gtf + ".gz"

            if os.path.exists(input_gtf) and os.path.exists(output_gtf):
                return output_gtf

            try:
                w = open(input_gtf, "w+")

                with open(old_input_gtf) as r:
                    for line in r:
                        if line.startswith("#"):
                            w.write(line)
                            continue

                        lines = line.split()

                        if len(lines) < 1:
                            continue

                        data.append(
                            GenomicLoci(
                                chromosome=lines[0],
                                start=lines[3],
                                end=lines[4],
                                strand=lines[6],
                                gtf_line=line
                            )
                        )

                for i in sorted(data, key=lambda x: [x.chromosome, x.start]):
                    w.write(i.gtf_line)

                w.close()
            except IOError as err:
                logger.error(f"could not sort gtf, because {err}")
                exit(err)

        return output_gtf

    def __load__(self, region: GenomicLoci):
        u"""
        Load transcripts inside of region
        :param region: target region
        :return:
        """

        if not region.chromosome.startswith("chr"):
            logger.info("Guess need 'chr'")
            region.chromosome = "chr" + region.chromosome
        else:
            logger.info("Guess 'chr' is redundant")
            region.chromosome = region.chromosome.replace("chr", "")

        transcripts = {}

        r = pysam.TabixFile(self.path)

        try:
            iter_ = r.fetch(region.chromosome, region.start, region.end, parser=pysam.asGTF())
        except ValueError:
            try:
                if not region.chromosome.startswith("chr"):
                    logger.info("Guess need 'chr'")
                    iter_ = r.fetch(
                        "chr" + region.chromosome,
                        region.start, region.end, parser=pysam.asGTF()
                    )
                else:
                    logger.info("Guess 'chr' is redundant")
                    iter_ = r.fetch(
                        region.chromosome.replace("chr", ""),
                        region.start, region.end, parser=pysam.asGTF()
                    )
            except ValueError as err:
                logger.warn("please check the input region and gtf files")
                logger.error(err)
                raise err

        for rec in iter_:
            if re.search(r"(rna|transcript|cds)", rec.feature, re.I):
                if rec.transcript_id not in transcripts.keys():
                    transcripts[rec.transcript_id] = Transcript(
                        chromosome=rec.contig,
                        start=max(rec.start + 1, region.start),
                        end=min(rec.end + 1, region.end),
                        strand=rec.strand,
                        transcript_id=rec.transcript_id,
                        gene_id=rec.gene_id,
                        gene=rec.gene_name,
                        transcript=rec.transcript_name,
                        exons=[]
                    )

            elif re.search(r"(exon)", rec.feature, re.I):
                if rec.start + 1 >= region.end or rec.end + 1 <= region.start:
                    continue

                transcripts[rec.transcript_id].exons.append(
                    GenomicLoci(
                        chromosome=rec.contig,
                        start=max(rec.start + 1, region.start),
                        end=min(rec.end + 1, region.end),
                        strand=rec.strand
                    )
                )

        self.transcripts = sorted(transcripts.values())
        r.close()


if __name__ == '__main__':
    pass

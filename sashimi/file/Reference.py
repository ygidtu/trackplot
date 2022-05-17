#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2020.05.07

This scripts contains the class handle the reference file
"""
import gzip
import os
import re
from typing import List, Optional

import filetype
import matplotlib as mpl
import pysam
from matplotlib import pyplot as plt

from conf.logger import logger
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.Readder import Reader
from sashimi.base.Transcript import Transcript
from sashimi.base.Protein import CdsProtein
from sashimi.file.File import File

# Put here to avoid outline stroke of font for this moment
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.family"] = 'Arial'


class Reference(File):
    u"""
    The reference file, support gtf and gff format
    """

    def __init__(self, path: str, category: str = "gtf"):
        u"""
        init func
        :param path: path to input file
        """
        categories = ["gtf", "bam"]
        assert category in categories, f"category should be one of {categories}, instead of {category}"

        super().__init__(path=self.index_gtf(path) if category == "gtf" else path)
        self.category = category
        self.data = []
        self.domain = None

    def __add__(self, other):
        assert isinstance(other, Reference), "only Reference and Reference could be added"
        new_ref = Reference(self.path, category=self.category)
        new_ref.data += self.data
        new_ref.data += other.data
        new_ref.data = sorted(new_ref.data)
        if self.domain:
            new_ref.__add_domain__()
        return new_ref

    @classmethod
    def create(cls, path: str, category: str = "gtf"):
        u"""
        create reference file object
        :param path: path to input file
        :param category: the type of reference file, include gtf and bam: customized reads as references
        :return: Reference obj
        """
        assert os.path.exists(path), f"{path} not exists"
        return cls(path=path, category=category)

    def __add_domain__(self):
        gene_id = set(map(lambda x: x.gene_id, self.data))
        transcript_id = set(map(lambda x: x.transcript_id, self.data))
        chromosome_id = self.data[0].chromosome

        # if domain is not None, then add it.
        if self.domain:
            gene_id = gene_id.difference(self.domain.gene_id)
            transcript_id = transcript_id.difference(self.domain.transcript_id)

        if len(transcript_id) != 0 and len(gene_id) != 0:
            with pysam.Tabixfile(self.path) as gtf_tabix:
                domain_info = CdsProtein.__re_iter_gtf__(
                    gtf_tabix=gtf_tabix,
                    chromosome=chromosome_id,
                    transcript_id=transcript_id,
                    gene_id=chromosome_id
                )
            # add pep information to cdsProtein.
            if self.domain:
                self.domain.add(domain_info)
            else:
                self.domain = domain_info

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
    def sort_gtf(cls, input_gtf: str, output_gtf: str):
        data = []
        logger.info("Sorting %s" % input_gtf)

        try:
            w = open(output_gtf.replace(".gz", ""), "w+")

            with open(input_gtf) as r:
                for line in r:
                    if line.startswith("#"):
                        w.write(line)
                        continue

                    lines = line.split()

                    if len(lines) < 1:
                        continue

                    data.append(
                        GenomicLoci(
                            chromosome=lines[0], start=lines[3], end=lines[4],
                            strand=lines[6], gtf_line=line
                        )
                    )

            for i in sorted(data, key=lambda x: [x.chromosome, x.start]):
                w.write(i.gtf_line)
            w.close()

            pysam.tabix_compress(output_gtf.replace(".gz", ""), output_gtf)
            os.remove(output_gtf.replace(".gz", ""))
        except IOError as err:
            logger.error(f"could not sort gtf, because {err}")
            exit(err)

    @classmethod
    def index_gtf(cls, input_gtf):
        u"""
        Check the tabix index of input gtf file

        Extract only exon tags and keep it clean
        :param input_gtf: path to input gtf file
        :return path to compressed and indexed bgzipped gtf file
        """
        gtf = cls.is_gtf(input_gtf)
        assert gtf % 10 == 1, f"{input_gtf} seems not be gtf format"

        index = False
        if gtf // 10 > 0:
            output_gtf = input_gtf
        else:
            output_gtf = input_gtf + ".gz"

        sorted_gtf = re.sub(r"\.gtf(.gz)?$", "", input_gtf) + ".sorted.gtf.gz"
        if os.path.exists(sorted_gtf) and os.path.exists(sorted_gtf + ".tbi"):
            return sorted_gtf

        if not os.path.exists(output_gtf) or not os.path.exists(output_gtf + ".tbi"):
            logger.info(f"the {output_gtf} or tbi index not exists")
            index = True
        elif os.path.getctime(output_gtf) < os.path.getctime(output_gtf):
            logger.info("the tbi index is older than the gtf file")
            index = True

        if not index:
            return output_gtf

        logger.info("Create index for %s", input_gtf)
        try:
            cls.sort_gtf(input_gtf, sorted_gtf)
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

            cls.sort_gtf(input_gtf, sorted_gtf)
            pysam.tabix_index(
                sorted_gtf, preset="gff",
                force=True, keep_original=True
            )

        return output_gtf

    def __load_gtf__(self, region: GenomicLoci) -> List[Transcript]:

        u"""
        Load transcripts inside of region from gtf file
        :param region: target region
        :return: list of Transcript
        """
        transcripts = {}

        for rec in Reader.read_gtf(self.path, region):
            start = max(rec.start, region.start)
            end = min(rec.end, region.end)

            if end + 1 <= region.start:
                continue
            if start + 1 >= region.end:
                break

            if re.search(r"(rna|transcript|cds)", rec.feature, re.I):
                if rec.transcript_id not in transcripts.keys():
                    transcripts[rec.transcript_id] = Transcript(
                        chromosome=rec.contig,
                        start=start,
                        end=end,
                        strand=rec.strand,
                        transcript_id=rec.transcript_id,
                        gene_id=rec.gene_id,
                        gene=rec.gene_name,
                        transcript=rec.transcript_name,
                        exons=[]
                    )

            elif re.search(r"(exon)", rec.feature, re.I):
                transcripts[rec.transcript_id].exons.append(
                    GenomicLoci(
                        chromosome=rec.contig,
                        start=start,
                        end=end,
                        strand=rec.strand
                    )
                )
        return sorted(transcripts.values())

    def __load_bam__(self, region: GenomicLoci, threshold_of_reads: int = 0) -> List[Transcript]:
        u"""
        Load transcripts inside of region from bam file
        :param region: target region
        :param threshold_of_reads: only kept reads with minimum frequency
        :return: list of Transcript
        """
        transcripts = {}
        try:
            for read, strand in Reader.read_bam(self.path, region):
                start = read.reference_start

                exons_in_read = []
                for cigar, length in read.cigartuples:
                    cur_start = start + 1
                    cur_end = start + length + 1

                    if cigar == 0:  # M
                        if cur_start < region.end < cur_end:
                            exons_in_read.append(GenomicLoci(
                                chromosome=read.reference_name,
                                start=cur_start if cur_start > region.start else region.start,
                                end=cur_end if cur_end <= region.end else region.end,
                                strand="+",
                            ))

                    elif cigar not in (1, 2, 4, 5):  # I, D, S, H
                        start += length

                t = Transcript(
                    chromosome=read.reference_name,
                    start=read.reference_start + 1 if read.reference_start + 1 > region.start else region.start,
                    end=read.reference_end + 1 if read.reference_end + 1 < region.end else region.end,
                    strand=strand,
                    exons=exons_in_read
                )

                # t.gene = str(t)
                # t.transcript = str(t)
                # t.gene_id = str(t)
                # t.transcript_id = str(t)
                transcripts[t] = transcripts.get(t, 0) + 1
        except IOError as err:
            logger.error('There is no .bam file at {0}'.format(self.path))
            logger.error(err)
        except ValueError as err:
            logger.error(self.path)
            logger.error(err)

        return sorted([x for x, y in transcripts.items() if y > threshold_of_reads])

    def load(self, region: GenomicLoci, domain: Optional[bool] = False, threshold_of_reads: int = 0):
        u"""
        Load transcripts inside of region
        :param region: target region
        :param domain: init with domain, default: False
        :param threshold_of_reads: used for bam file, only kept reads with minimum frequency
        :return:
        """
        assert isinstance(region, GenomicLoci), "region should be a GenomicLoci object"
        self.region = region
        if self.category == "gtf":
            self.data = self.__load_gtf__(region)

            if domain:
                self.__add_domain__()
        else:
            self.data = self.__load_bam__(region, threshold_of_reads)

    def add_interval(self, interval_file: str, interval_label: str):
        u"""
        Add another annotation in to data, and each annotation has a track
        :param interval_file: a bed file
        :return:
        """
        try:
            if not os.path.exists(interval_file + ".tbi"):
                interval_file = pysam.tabix_index(
                    interval_file,
                    preset="bed",
                    force=True,
                    keep_original=True
                )

            interval_target = []
            for rec in Reader.read_gtf(interval_file, region=self.region, bed=True):
                start = max(rec.start, self.region.start)
                end = min(rec.end, self.region.end)

                if end + 1 <= self.region.start:
                    continue
                if start + 1 >= self.region.end:
                    break

                try:
                    strand = rec.strand
                    rec_name = rec.name
                except KeyError:
                    strand = "*"
                    rec_name = ""

                interval_target.append(
                    GenomicLoci(
                        chromosome=rec.contig,
                        start=rec.start,
                        end=rec.end,
                        strand=strand,
                        name=rec_name
                    )
                )

            if len(interval_target) != 0:
                self.data.append(Transcript(
                    chromosome=rec.contig,
                    start=start,
                    end=end,
                    strand=strand,
                    exons=interval_target,
                    transcript_id=interval_label,
                    gene_id=interval_label,
                    gene=interval_label,
                    transcript=interval_label,
                    category="interval"
                ))

        except NameError:
            raise "Target region was not found, run load before add_bed."

        except OSError:
            raise f"Error found when build index for {interval_file}, please sort file manually"


if __name__ == "__main__":
    # region = GenomicLoci("chr1", 1270656, 1284730, "+")
    # print(len(region))
    # ref = Reference.create("../../example/example.gtf")
    # ref.load(region)
    # print(len(ref.data))
    #
    # ref1 = Reference.create("../../example/bams/1.bam", category="bam")
    # ref1.load(region, 10)
    # print(len(ref1.data))

    loc = GenomicLoci(chromosome="chr1",
                      start=1017198,
                      end=1051741,
                      strand="-")

    gtf_ref = Reference("../../example/example.sorted.gtf.gz")
    gtf_ref.load(loc, domain=True)
    print(gtf_ref.domain)
    # gene_id = set(map(lambda x: x.gene_id, gtf_ref.transcripts))
    # transcript_id = set(map(lambda x: x.transcript_id, gtf_ref.transcripts))
    # print(gene_id)
    # print(transcript_id)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import gzip

import random

from cmdbench import benchmark_command, BenchmarkResults
from subprocess import check_call, check_output

import pysam
from tqdm import tqdm


class Reference(object):

    def __init__(self, line: str, id: str = None):
        self.line = line
        self.id = id
        rec = line.strip().split("\t")

        self.contig = rec[0]
        self.start = int(rec[3])
        self.end = int(rec[4])

    def __str__(self):
        return self.line

    def __hash__(self):
        return hash(self.line)

    def __lt__(self, other):
        if self.contig != other.contig:
            return self.contig < other.contig
        if self.start != other.start:
            return self.start < other.start
        return self.end < other.end

    def __gt__(self, other):
        if self.contig != other.contig:
            return self.contig > other.contig
        if self.start != other.start:
            return self.start > other.start
        return self.end > other.end

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

class Bench(object):

    __slots__ = ["res", "debug", "n", "init_time", "time", "memory"]

    def __init__(self, debug: bool = False):
        self.res = BenchmarkResults()
        self.debug = debug
        self.n = 0
        self.init_time = 0
        self.time = 0
        self.memory = 0

    def set_init_time(self, bench):
        self.init_time = bench.time

    def add(self, cmd, with_activate: bool = False):
        # print(cmd)
        if self.debug:
            check_call(cmd, shell=True)
        else:
            self.res.add_benchmark_result(benchmark_command(f"bash -c '{cmd}'", iterations_num=1))

            if self.res.iterations[-1]["process"]["exit_code"] != 0:
                print(self.res.iterations[-1]["process"]["stdout_data"])
                print(self.res.iterations[-1]["process"]["stderr_data"])
                raise ValueError(f"exit_code != 0: {cmd}")

            self.time = self.time + self.res.iterations[-1]['process']['execution_time']
            self.memory = max(self.memory, self.res.iterations[-1]['memory']['max'])
        if with_activate:
            self.n += 1

    def stats(self):
        return self.time - self.init_time * self.n, self.memory

def get_env(env: str):
    out = check_output("conda env list", shell=True).decode("utf-8")
    envs = {}
    for line in out.split("\n"):
        line = line.split()
        if len(line) == 2:
            envs[line[0]] = line[1]
    return envs.get(env, "")


def format_attributes(rec: str):
    attrs = {}
    for attr in rec.split(";"):
        attr = attr.split() if "=" not in attr else attr.split("=")
        if len(attr) > 1:
            attr = [x.strip().replace("'", "").replace('"', '') for x in attr]
            attrs[attr[0]] = attr[-1]
    return attrs


def subset_gtf_or_gff3(path, output, is_gtf=True, gene_ids=10000, seed=1):
    genes, transcripts, exons = {}, {}, {}

    r = gzip.open(path, "rt") if path.endswith(".gz") else open(path)

    for line in tqdm(r, desc="Reading..."):
        if line.startswith("#"):
            continue

        line = line.strip()
        lines = line.split("\t")
        if len(lines) < 8:
            continue
        attrs = format_attributes(lines[8])

        if is_gtf:
            gene_id, transcript_id, exon_id = attrs.get("gene_id"), attrs.get("transcript_id"), attrs.get("exon_number")

            # gene
            if gene_id and gene_id not in genes.keys():
                genes[gene_id] = Reference(line, gene_id)

            # transcript
            elif transcript_id and transcript_id not in transcripts.keys():
                if gene_id not in transcripts:
                    transcripts[gene_id] = []
                transcripts[gene_id].append(Reference(line, transcript_id))
            elif transcript_id:
                if transcript_id not in exons:
                    exons[transcript_id] = []
                exons[transcript_id].append(Reference(line, exon_id))
        else:
            if "Parent" not in attrs.keys() and lines[2] == "gene":
                gene_id = attrs.get("ID")
                if gene_id:
                    genes[gene_id] = Reference(line, gene_id)
            else:
                parent_id = attrs.get("Parent")

                if parent_id in genes.keys():
                    gene_id = parent_id
                    transcript_id = attrs.get("ID")

                    if gene_id not in transcripts:
                        transcripts[gene_id] = []
                    transcripts[gene_id].append(Reference(line, transcript_id))
                else:
                    exon_id = attrs.get("ID")
                    transcript_id = attrs.get("Parent")

                    if transcript_id not in exons:
                        exons[transcript_id] = []
                    exons[transcript_id].append(Reference(line, exon_id))

    r.close()

    # filter genes
    genes_for_selected = []
    for g in genes.keys():
        if len(transcripts.get(g, [])) > 2:
            genes_for_selected.append(g)
    genes_for_selected = sorted(genes_for_selected)

    if isinstance(gene_ids, int):
        if len(genes) > gene_ids:
            random.seed(seed)
            gene_ids = random.choices(genes_for_selected, k=gene_ids)
        else:
            gene_ids = genes_for_selected

    genes = [genes[x] for x in gene_ids]

    res = []

    with open(output, "w+") as w:
        for gene in tqdm(sorted(genes), desc="Saving..."):
            if is_gtf:
                res.append(gene)
            else:
                w.write(f"{gene}\n")
            for transcript in sorted(transcripts.get(gene.id, [])):
                if is_gtf:
                    res.append(transcript)
                else:
                    w.write(f"{transcript}\n")
                for exon in sorted(exons.get(transcript.id, [])):
                    if is_gtf:
                        res.append(exon)
                    else:
                        w.write(f"{exon}\n")

        for r in sorted(res):
            w.write(f"{r}\n")

    if is_gtf:
        pysam.tabix_compress(output, output + ".gz", force=True)
        pysam.tabix_index(output + ".gz", force=True, preset="gff")

    return gene_ids


def generate_event_region(path: str, key: str = "gene:ENSG00000186007"):
    with pysam.Tabixfile(path) as r:
        for rec in r.fetch(parser=pysam.asTuple()):
            if len(rec) < 8 or "gene" not in rec[2]:
                continue
            contig, start, end, strand, attributes = rec[0], int(rec[3]), int(rec[4]), rec[6], rec[8]

            if key in attributes:
                return contig, start - 1, end + 1, strand

    raise ValueError(f"There is no {key} in {path}")


def generate_gff3(path: str, output: str, key: str = "gene:ENSG00000186007"):
    contig, start, end, _ = generate_event_region(path, key)
    with pysam.Tabixfile(path) as r:
        if contig and start and end:
            with open(output, "w+") as w:
                for rec in r.fetch(contig, start - 1, end + 1):
                    w.write(rec + "\n")
        else:
            raise ValueError(f"There is no {key} in {path}")


if __name__ == '__main__':
    pass

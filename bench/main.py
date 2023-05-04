#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

import click
import gzip
import pysam

from run_ggsashimi import run_ggsashimi
from run_miso import run_miso, run_miso_preprocess
from run_trackplot import run_trackplot
from bench import subset_gtf_or_gff3


class ProjectStruct(object):

    def __init__(self, path: str):
        os.makedirs(path, exist_ok=True)
        self.stats = os.path.join(path, "stats.txt")

        self.funcs = {
            "trackplot": run_trackplot,
            "ggsashimi": run_ggsashimi,
            "miso": run_miso,
        }

        self.path = {x: os.path.join(path, x) for x in self.funcs.keys()}

        for i in self.path.values():
            os.makedirs(i, exist_ok=True)

        self.bam_lists = {k: os.path.join(v, "bam.list") for k, v in self.path.items()}

    def __iter__(self):
        for k, v in self.funcs.items():
            yield k, v

    def output(self, key: str, event: str):
        return os.path.join(self.path[key], event)


class FileList(object):

    def __init__(self, line: str):
        line = line.strip().split()
        headers = ["path", "key"]

        self.data = {x: y for x, y in zip(headers, line)}
        self.data["additional"] = line[len(headers):]

    @property
    def key(self):
        return self.data["key"]

    @property
    def path(self):
        assert os.path.exists(self.data['path']), f"{self.data['path']} not exists."
        return os.path.abspath(self.data["path"])

    @property
    def color(self):
        return self.data.get("color")

    @property
    def additional(self):
        if len(self.data["additional"]) > 0:
            return "\t".join(self.data["additional"])
        return ""

    def to_str(self, ID: int, format_: str = "trackplot"):
        res = ""
        if format_ == "ggsashimi":
            res = f"{self.key}_{ID}\t{self.path}"
        elif format_ == "trackplot":
            res = f"{self.path}\tbam\t{self.key}_{ID}"
        elif format_ == "miso":
            res = f"{self.path}\t{self.key}_{ID}"

        if self.additional:
            res = f"{res}\t{self.additional}"
        return res


@click.command(context_settings=dict(help_option_names=['-h', '--help']), no_args_is_help=True)
@click.option("-i", "--infile", type=click.Path(exists=True),
              help="Path to filelist with three columns, 1st path to bam; 2nd alias of bam; 3rd color of bam.",
              required=True)
@click.option("-o", "--output", type=click.Path(), help="Path to output directory.", required=True)
@click.option("-g", "--reference", type=str, help="Prefix of gtf and gff3.")
@click.option("-e", "--event", type=int, help="The number of compared genes", default=10000, show_default=True)
@click.option("-r", "--repeat", type=click.IntRange(min=1), help="How many files to generated.",
              default=6, show_default=True)
@click.option("-n", "--n-jobs", type=click.IntRange(min=1), help="How many processes to use.",
              default=1, show_default=True)
@click.option("-a", "--append", is_flag=True, help="Append new results to exist file.",
              show_default=True)
@click.option("-s", "--seed", type=int, default=1, help="The seed for random selection of genes.", show_default=True)
@click.option("-s", "--seed", type=int, default=1, help="The seed for random selection of genes.", show_default=True)
def main(infile, output, repeat, n_jobs, reference, event, seed, append):
    data = []
    with open(infile) as r:
        for line in r:
            data.append(FileList(line))

    proj = ProjectStruct(output)

    if not append:
        for key, path in proj.bam_lists.items():
            print(f"generating bam list of {key}")
            with open(path, "w+") as w:
                for i in range(repeat):
                    w.write(data[i % len(data)].to_str(ID=i, format_=key) + "\n")

        if not os.path.exists(f"{reference}.gff3") and os.path.exists(f"{reference}.gff3.gz"):
            print(f"generating {reference}.gff3")
            with open(f"{reference}.gff3", "w+") as w:
                with gzip.open(f"{reference}.gff3.gz", "rt") as r:
                    for line in r:
                        w.write(line)

        print("down sampling gff3")
        gene_ids = subset_gtf_or_gff3(
            f"{reference}.gff3", os.path.join(output, "target.gff3"),
            is_gtf=False, gene_ids=event, seed=seed
        )

        print("down sampling gtf")
        subset_gtf_or_gff3(
            f"{reference}.gtf", os.path.join(output, "target.gtf"),
            is_gtf=True, gene_ids={x.split(":")[-1] for x in gene_ids}, seed=seed
        )

    events = {}
    fcode = "a+" if append else "w+"

    finished = set()
    if append:
        with open(proj.stats, "r") as r:
            for line in r:
                line = line.split()
                finished.add(line[0])

    with open(proj.stats, fcode) as w:
        if not append:
            w.write("event\ttime\tmemory\tsoftware\tnum_of_files\tn_jobs\n")

        p_t, p_m = run_miso_preprocess(**{
            "output": proj.path["miso"],
            "gff": os.path.join(output, "target.gff3"),
            "bam": proj.bam_lists["miso"],
            "env": "misopy",
            "n_jobs": n_jobs,
        })
        print(f"preprocess\t{p_t}\t{p_m}\tmiso_preprocess\t{repeat}\t{n_jobs}")
        w.write(f"preprocess\t{p_t}\t{p_m}\tmiso_preprocess\t{repeat}\t{n_jobs}\n")

        gene_ids_pkl = []
        for parent, _, files in os.walk(proj.path["miso"]):
            for f in files:
                if f.endswith(".pickle"):
                    if append and f.replace(".pickle", "") in finished:
                        continue
                    gene_ids_pkl.append(f.replace(".pickle", ""))

        for e in sorted(set(gene_ids_pkl)):
            e_id = 0
            while f"{e}.{e_id}" in events.keys():
                e_id += 1

            for key, func in proj:
                try:
                    t, m = func(**{
                        "event": e if "miso" == key else e.split(":")[-1],
                        "output": proj.output(key, "plots") if "miso" != key else proj.path[key],
                        "gtf": os.path.join(output, "target.gtf.gz"),
                        "bam": proj.bam_lists[key],
                        "env": key if "miso" != key else "misopy",
                        "n_jobs": n_jobs,
                    })
                    print(f"{e}\t{t}\t{m}\t{key}\t{repeat}\t{n_jobs}")
                    w.write(f"{e}\t{t}\t{m}\t{key}\t{repeat}\t{n_jobs}\n")
                    w.flush()
                except ValueError as err:
                    print(err)


if __name__ == '__main__':
    # python main.py -i bam_list.txt -o benchmark_files/ -a -g ref/Homo_sapiens.GRCh38.101.chr --repeat 1 --n-jobs 1 --event 3
    main()

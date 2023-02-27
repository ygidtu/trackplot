#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

import click
import gzip

from run_ggsashimi import run_ggsashimi
from run_miso import run_miso
from run_trackplot import run_trackplot


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
@click.option("-e", "--event", type=str, help="comma seperated gene ids",
              default="ENSG0000022397,ENSG0000022723", show_default=True)
@click.option("-r", "--repeat", type=click.IntRange(min=1), help="How many files to generated.",
              default=6, show_default=True)
@click.option("-n", "--n-jobs", type=click.IntRange(min=1), help="How many processes to use.",
              default=1, show_default=True)
@click.option("-a", "--append", is_flag=True, help="Append new results to exist file.",
              show_default=True)
def main(infile: str, output: str, repeat: int, n_jobs: int, reference: str, event: str, append: bool):
    for postfix in [".gff3.gz", ".gff3.gz.tbi", ".gtf.gz", ".gtf.gz.tbi"]:
        assert os.path.exists(f"{reference}{postfix}"), f"{reference}{postfix} not exists"

    data = []
    with open(infile) as r:
        for line in r:
            data.append(FileList(line))

    proj = ProjectStruct(output)
    for key, path in proj.bam_lists.items():
        print(f"generating bam list of {key}")
        with open(path, "w+") as w:
            for i in range(repeat):
                w.write(data[i % len(data)].to_str(ID=i, format_=key) + "\n")

    events = {}

    if not os.path.exists(f"{reference}.gff3") and os.path.exists(f"{reference}.gff3.gz"):
        with open(f"{reference}.gff3", "w+") as w:
            with gzip.open(f"{reference}.gff3.gz", "rt") as r:
                for line in r:
                    w.write(line)

    fcode = "a+" if append else "w+"
    with open(proj.stats, fcode) as w:
        if not append:
            w.write("event\ttime\tmemory\tsoftware\tnum_of_files\tn_jobs\n")
        for e in event.split(","):
            e_id = 0
            while f"{e}.{e_id}" in events.keys():
                e_id += 1

            for key, func in proj:
                t, m = func(**{
                    "event": e if "miso" != key else f"gene:{e}",
                    "output": proj.output(key, f"{e}.{e_id}"),
                    "gff": f"{reference}.gff3.gz", # if key == "miso" else f"{reference}.gff3",
                    "gtf": f"{reference}.gtf.gz",
                    "bam": proj.bam_lists[key],
                    "env": key if "miso" != key else "misopy",
                    "n_jobs": n_jobs, "generate_gff": key == "miso"
                })
                print(f"{e}\t{t}\t{m}\t{key}\t{repeat}\t{n_jobs}")
                w.write(f"{e}\t{t}\t{m}\t{key}\t{repeat}\t{n_jobs}\n")
                w.flush()


if __name__ == '__main__':
    main()

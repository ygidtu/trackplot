#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

import click

from run_ggsashimi import run_ggsashimi
from run_miso import run_miso
from run_sashimipy import run_sashimipy


class ProjectStruct(object):

    def __init__(self, path: str):
        os.makedirs(path, exist_ok=True)
        self.stats = os.path.join(path, "stats.txt")

        self.path = {
            "sashimipy": os.path.join(path, "sashimipy"),
            "ggsashimi": os.path.join(path, "ggsashimi"),
            "miso": os.path.join(path, "miso")
        }

        for i in self.path.values():
            os.makedirs(i, exist_ok=True)

        self.bam_lists = {k: os.path.join(v, "bam.list") for k, v in self.path.items()}

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

    def to_str(self, ID: int, format_: str = "sashimipy"):
        res = ""
        if format_ == "ggsashimi":
            res = f"{self.key}_{ID}\t{self.path}"
        elif format_ == "sashimipy":
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
def main(infile: str, output: str, repeat: int, n_jobs: int, reference: str, event: str):
    for postfix in [".gff3.gz", ".gff3.gz.tbi", ".gtf.gz", ".gtf.gz.tbi"]:
        assert os.path.exists(f"{reference}{postfix}"), f"{reference}{postfix} not exists"

    # if os.path.exists(path):
    #     rmtree(path)

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

    with open(proj.stats, "w+") as w:
        w.write("event\ttime\tmemory\tsoftware\tnum_of_files\tn_jobs\n")
        for e in event.split(","):

            for key, func in {
                "miso": run_miso,
                "sashimipy": run_sashimipy,
                "ggsashimi": run_ggsashimi
            }.items():
                t, m = func(**{
                    "event": e if key != "miso" else f"gene:{e}",
                    "output": proj.output(key, e),
                    "gff": f"{reference}.gff3.gz",
                    "gtf": f"{reference}.gtf.gz",
                    "bam": proj.bam_lists[key],
                    "env": key if key != "miso" else "misopy",
                    "n_jobs": n_jobs
                })

                w.write(f"{e}\t{t}\t{m}\t{key}\t{repeat}\t{n_jobs}\n")


if __name__ == '__main__':
    main()

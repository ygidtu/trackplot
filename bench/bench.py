#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from subprocess import check_call, check_output

import pysam
from cmdbench import benchmark_command, BenchmarkResults


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

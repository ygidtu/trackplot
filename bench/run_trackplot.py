#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""

"""
import os
from shutil import rmtree

import click

from bench import Bench, generate_event_region, get_env


def run_trackplot(event: str, output: str, gtf: str, bam: str, env: str, n_jobs: int, **kwargs):
    bench = Bench()

    root = get_env(env)
    trackplot_plot = "trackplot"
    with_activate = False
    if root:
        print(f"Using conda env: {root}")
        trackplot_plot = f"source activate {env} && {trackplot_plot}"
        activate = Bench()
        activate.add(f"source activate {env}")
        bench.set_init_time(activate)
        with_activate = True

    contig, start, end, strand = generate_event_region(gtf, event)
    os.makedirs(output, exist_ok=True)

    with open(os.path.join(output, "bam.list"), "w+") as w:
        with open(bam) as r:
            for line in r:
                line = line.strip().split("\t")
                assert os.path.exists(line[0]), f"{line[0]} not exists"
                line[0] = os.path.abspath(line[0])
                w.write("\t".join(line) + "\n")

    bench.add(f"{trackplot_plot} --event {contig}:{start}-{end}:{strand} "
              f"--density {os.path.join(output, 'bam.list')} -p {n_jobs} -r {gtf} "
              f"--output {os.path.join(output, event)}.pdf --height 0.01", with_activate=with_activate)
    stat = bench.stats()
    del bench
    return stat


@click.command(context_settings=dict(help_option_names=['-h', '--help']), no_args_is_help=True)
@click.option("-e", "--event", type=str, help="The target gene ID.", required=True)
@click.option("-g", "--gtf", type=click.Path(exists=True), help="Path to gtf.", required=True)
@click.option("-b", "--bam", type=click.Path(exists=True),
              help="List of bams for trackplot.", required=True)
@click.option("-o", "--output", type=click.Path(), help="Path to output dir.", required=True)
@click.option("--env", type=str, help="Name of used conda env.")
@click.option("--n-jobs", type=int, default=1, help="The number of processes to use.", show_default=True)
def main(**kwargs):
    if os.path.exists(kwargs["output"]):
        rmtree(kwargs["output"])
    execution_time, max_memory = run_trackplot(**kwargs)
    print(f"execution_time={execution_time}; max_memory={max_memory}")


if __name__ == '__main__':
    main()

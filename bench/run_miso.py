#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""

"""
import configparser
import os
from shutil import rmtree
import click
import pysam

from bench import Bench, generate_gff3, get_env


class Config(object):
    def __init__(self, bam_prefix: str, miso_prefix: str):
        self.config = configparser.ConfigParser()
        self.bam_prefix = bam_prefix
        self.miso_prefix = miso_prefix
        self.colors = []
        self.bams = []
        self.miso = []

    def write(self, path: str):
        self.config['data'] = {
            "bam_prefix": self.bam_prefix,
            "miso_prefix": self.miso_prefix,
            "bam_files": self.bams,
            "miso_files": self.miso
        }

        self.config['plotting'] = __PLOTTING__ = {
            'fig_width': 7, 'fig_height': 5,
            'intron_scale': 30, 'exon_scale': 4,
            'logged': False, 'font_size': 6, 'ymax': 150,
            'show_posteriors': True, 'bar_posteriors': False,
            'number_junctions': True, 'resolution': .5,
            'posterior_bins': 40, 'gene_posterior_ratio': 5,
            'colors': self.colors
        }

        with open(path, "w+") as w:
            self.config.write(w)

    def add_file(self, bam_name: str, miso_name: str, colors: str):
        self.bams.append(bam_name)
        self.miso.append(miso_name)
        self.colors.append(colors)


class Files:

    def __init__(self, output):
        output = os.path.abspath(output)
        os.makedirs(output, exist_ok=True)
        self.gtf = os.path.join(output, "target.gff3")
        self.index = os.path.join(output, "index")
        self.miso = os.path.join(output, "miso")
        self.conf = os.path.join(output, "plot_config.ini")
        self.plot = os.path.join(output, "plots")


def run_miso(event: str, output: str, gff: str, bam: str, env: str, n_jobs: int, **kwargs):
    bench = Bench()

    root = get_env(env)
    index_gff = "index_gff"
    miso = "miso"
    sashimi_plot = "sashimi_plot"
    if root:
        print(f"Using conda env: {root}")
        index_gff = f"source activate {env} && index_gff"
        miso = f"source activate {env} && miso"
        sashimi_plot = f"source activate {env} && sashimi_plot"
        activate = Bench()
        activate.add(f"source activate {env}")
        bench.set_init_time(activate)

    f = Files(output)
    generate_gff3(gff, f.gtf, event)
    bench.add(f"{index_gff} --index {f.gtf} {f.index}", with_activate=True)

    config = None
    bams = {}
    with open(bam) as r:
        for line in r:
            path, key, color = line.strip().split()
            assert os.path.exists(path), f"{path} not exists"
            if config is None:
                config = Config(bam_prefix=os.path.dirname(os.path.abspath(path)), miso_prefix=f.miso)

            config.add_file(os.path.basename(path), miso_name=key, colors=color)
            bams[key] = path
            o = os.path.join(f.miso, key)
            os.makedirs(o, exist_ok=True)

            read_len = 44
            with pysam.AlignmentFile(path) as r:
                for rec in r:
                    read_len = rec.infer_read_length()
                    break

            bench.add(f"{miso} --run {f.index} {path} --output-dir {o} --read-len {read_len} -p {n_jobs}", with_activate=True)

    config.write(f.conf)

    os.makedirs(f.plot, exist_ok=True)
    bench.add(f"{sashimi_plot} --plot-event {event} {f.index} {f.conf} --output-dir {f.plot}")
    return bench.stats()


@click.command(context_settings=dict(help_option_names=['-h', '--help']), no_args_is_help=True)
@click.option("-e", "--event", type=str, help="The target gene ID.", required=True)
@click.option("-g", "--gff", type=click.Path(exists=True), help="Path to gff3.", required=True)
@click.option("-b", "--bam", type=click.Path(exists=True),
              help="Path to bam list, tab seperated list, 1st column is path to bam, "
                   "2nd column is alias of bam, 3rd is color.", required=True)
@click.option("-o", "--output", type=click.Path(), help="Path to output dir.", required=True)
@click.option("--env", type=str, help="Name of used conda env.", default="misopy")
@click.option("--n-jobs", type=int, default=1, help="The number of processes to use.", show_default=True)
def main(**kwargs):
    if os.path.exists(kwargs["output"]):
        rmtree(kwargs["output"])
    execution_time, max_memory = run_miso(**kwargs)
    print(f"execution_time={execution_time}; max_memory={max_memory}")


if __name__ == '__main__':
    main()

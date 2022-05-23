#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""
from typing import List, Optional, Set, Union, Dict
from copy import deepcopy

import matplotlib.pyplot as plt

from matplotlib import gridspec

from conf.logger import logger
from sashimi.base.Stroke import Stroke
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.file.Reference import Reference
from sashimi.file.File import File
from sashimi.file.Fasta import Fasta
from sashimi.file.Bam import Bam
from sashimi.file.Bigwig import Bigwig
from sashimi.file.Depth import Depth
from sashimi.base.ReadDepth import ReadDepth
from sashimi.plot_func import plot_line, plot_density, plot_reference, plot_heatmap, init_graph_coords, set_x_ticks, set_indicator_lines, set_focus, plot_stroke


class PlotInfo(object):
    u"""
    this class is used to collect all the plot information.
    """

    def __init__(self, obj: File, category: str = "", type_: str = "", group: str = ""):
        u"""
        init this class
        :param obj: the input file object
        :param category: the input file category
        :param type_: the plot type
        :param group: group name of input file, only used to plots type_ == heatmap
        """
        self.obj = [obj]
        self.group = group
        self.type = type_
        self.category = [category]

    def __str__(self) -> str:
        return ";".join([obj.path for obj in self.obj])

    @property
    def data(self) -> Dict[str, ReadDepth]:
        data = {}
        for obj in self.obj:
            if isinstance(obj.data, dict):
                data.update(data)
            else:
                data[obj.label] = obj.data
        return data

    def len(self, show_side_plot: bool = False) -> int:
        n = 0
        if not self.category:
            pass
        elif self.type not in ["heatmap", "line"] and self.category[0] == "bam" and show_side_plot:
            n += 2
        else:
            n += 1
        return n

    def add(self, obj: File, category: str = "", type_: str = ""):
        u"""
        add new input file to specific group
        """
        assert type_ in ["heatmap", "line"], "only heatmap/line plot needs add"

        self.obj.append(obj)
        self.category.append(category)
        return self

    def load(self, region:GenomicLoci, *args, **kwargs):
        for obj in self.obj:
            obj.load(region=region, *args, **kwargs)
        return self


class Plot(object):
    u"""
    this class is the main framework of sashimi
    """

    def __init__(self):
        u"""
        init this class
        """
        self.region = None
        self.sites = {}
        self.focus = {}
        self.stroke = []
        self.events = None
        self.sequence = None
        self.reference = None
        self.graph_coords = None
        self.plots = []

    @property
    def chrom(self) -> Optional[str]:
        if self.region:
            return self.region.chrom

    @property
    def start(self) -> Optional[int]:
        if self.region:
            return self.region.start

    @property
    def end(self) -> Optional[int]:
        if self.region:
            return self.region.end

    @property
    def strand(self) -> Optional[str]:
        if self.region:
            return self.region.strand

    @property
    def exons(self) -> Optional[List[List[int]]]:
        if self.reference:
            return self.reference.exons

    def set_region(self, chromosome: str, start: int, end: int, strand: str = "+"):
        u"""
        change the plot region
        :param chromosome:
        :param start:
        :param end:
        :param strand:
        :return:
        """
        self.region = GenomicLoci(chromosome, start=start, end=end, strand=strand)
        logger.info(f"set region to {self.region}")
        return self

    def add_sites(self, sites):
        u"""
        highlight specific sites
        :param sites: string in 100,200 format or int
        :return:
        """
        assert self.region is not None, f"please set plot region first."
        if sites is not None:
            if isinstance(sites, str):
                sites = [int(x) for x in sites.split(",")]

                for s in sites:
                    if s not in self.sites.keys():
                        self.sites[s - self.start] = "blue"
                    else:
                        self.sites[s - self.start] = "red"
            elif isinstance(sites, int):
                if sites not in self.sites.keys():
                    self.sites[sites - self.start] = "blue"
                else:
                    self.sites[sites - self.start] = "red"
        return self

    def add_focus(self, focus: Optional[str], start: int = 0, end: int = 0):
        u"""
        set focus region
        :param focus: string in 100-200:300-400
        :param start: start site
        :param end: end site
        :return:
        """
        assert self.region is not None, f"please set plot region first."
        if focus:
            for site in focus.split(":"):
                site = sorted([int(x) - self.start for x in site.split("-")])
                if site[0] < 0:
                    site[0] = 0

                if site[-1] > len(self):
                    site[-1] = len(self)

                self.focus[site[0]] = max(site[1], self.focus.get(site[0], -1))

        if 0 < start < end:
            self.focus[start] = max(end, self.focus.get(start, -1))
        return self

    def add_stroke(
            self,
            stroke: Optional[str] = None,
            start: int = 0,
            end: int = 0,
            label: str = "",
            color: str = "black"
    ):
        u"""
        add stroke to plot
        :param stroke: stroke string in 100-200@red;300-400 format
        :param start: start site of stroke
        :param end: end site of stroke
        :param label: label of stroke
        :param color: color of stroke
        :return:
        """
        assert self.region is not None, f"please set plot region first."
        if stroke:
            self.stroke += Stroke.create(stroke, self.region)

        if 0 < start < end:
            self.stroke.append(Stroke(start - self.start, end - self.end, color, label))
        return self

    def set_sequence(self, fasta: str):
        u"""
        set sequence info for
        :param fasta: path to indexed fasta file
        :return:
        """
        logger.info(f"fetch sequence from {fasta}")
        self.sequence = Fasta.create(fasta)
        return self

    def set_reference(self, gtf: str, interval: Optional[str] = None, interval_label: Optional[str] = None):
        u"""
        add transcripts to this region
        :param gtf:
        :param interval:
        :param interval_label:
        :return:
        """
        logger.info(f"set reference file to {gtf}")
        self.reference = Reference.create(gtf)

        if interval and interval_label:
            self.reference.add_interval(interval, interval_label)

        return self

    def add_interval(self, interval: str, interval_label: str):
        u"""

        """
        assert self.reference is not None, "please set_reference first."
        self.reference.add_interval(interval, interval_label)
        return self

    @staticmethod
    def __init_input_file__(path: str,
                            category: str = "bam",
                            label: Union[str, List[str]] = "",
                            title: str = "",
                            barcodes: Optional[Set[str]] = None,
                            cell_barcode: str = "BC",
                            umi_barcode: str = "UB",
                            library: str = "fr-unstrand"
                            ):
        if category == "bam":
            obj = Bam.create(
                path,
                label=label,
                title=title,
                barcodes=barcodes,
                cell_barcode=cell_barcode,
                umi_barcode=umi_barcode,
                library=library
            )
        elif category == "bigwig" or category == "bw":
            category = "bw"
            obj = Bigwig.create(path, label=label, title=title)
        elif category == "depth":
            obj = Depth.create(path, label=label, title=title)
        else:
            raise ValueError(f"the category should be one of [bam, bigwig, bw, depth], instead of {category}")
        return obj, category

    def add_density(self, path: str, category: str, *args, **kwargs):
        u"""
        add density object to plot
        :param path: the path to input file
        :param category: the input file type
        :return:
        """
        obj, category = self.__init_input_file__(path=path, category=category, *args, **kwargs)
        self.plots.append(PlotInfo(obj=obj, type_="density", category=category))
        return self

    def add_heatmap(self, path: str, group: str, category: str = "bam", *args, **kwargs):
        u"""
        add multiple objects for a group of heatmap
        :param path: path to input files
        :param group: the heatmap group
        :param category: file category corresponding to input file
        :return:
        """
        obj, category = self.__init_input_file__(path=path, category=category, *args, **kwargs)

        exists = False
        for p in self.plots:
            if p.group == group and p.type == "heatmap":
                p.add(obj=obj, category=category, type_="heatmap")
                exists = True
                break

        if not exists:
            self.plots.append(PlotInfo(obj=obj, category=category, type_="heatmap", group=group))

        return self

    def add_line(self, path: str, group: str, category: str = "bam", *args, **kwargs):
        u"""
        add multiple objects for a group of heatmap
        :param path: path to input files
        :param group: the heatmap group
        :param category: file category corresponding to input file
        :return:
        """
        obj, category = self.__init_input_file__(path=path, category=category, *args, **kwargs)

        exists = False
        for p in self.plots:
            if p.group == group and p.type == "line":
                p.add(obj=obj, category=category, type_="line")
                exists = True
                break

        if not exists:
            self.plots.append(PlotInfo(obj=obj, category=category, type_="line", group=group))

        return self

    def __len__(self) -> int:
        return self.end - self.start + 1

    def copy(self):
        return deepcopy(self)

    def plot(self,
             output: Optional[str] = None,
             show_side_plot: bool = False,
             reference_scale: Union[int, float] = .25,
             stroke_scale: Union[int, float] = .25,
             dpi: int = 300,
             fig_width: Union[int, float] = 0,
             fig_height: Union[int, float] = 0,
             *args, **kwargs):
        u"""
        draw image

        :param output: if output is empty then show this image by plt.showfig
        :param show_side_plot: whether to show side plot
        :param reference_scale: to adjust the size of reference plot
        :param stroke_scale: to adjust the size of stroke plot
        :param dpi: the dpi of saved plot
        :param fig_width: the width of figure, if width == 0, the let matplotlib decide the size of image
        :param fig_height: the height of figure, if height == 0, the let matplotlib decide the size of image
        """
        assert self.region is not None, f"please set the plotting region first."

        plots_n_rows = 0
        plots_n_cols = 1
        if self.reference is not None:
            logger.info("load reference")
            self.reference.load(self.region, *args, **kwargs)
            plots_n_rows += self.reference.len(scale=reference_scale)

        if self.stroke:
            plots_n_rows += int(max(len(self.stroke) * stroke_scale, 1))

        if self.sequence is not None:
            logger.info("load sequence")
            self.sequence.load(self.region)

        logger.info(f"load data")
        for p in self.plots:
            assert isinstance(p, PlotInfo), f"unrecognized data type: {type(p)}"
            try:
                p.load(region=self.region, *args, **kwargs)
            except Exception as err:
                logger.warning(f"failed to load data from {p}")
                raise err
            plots_n_rows += p.len(show_side_plot=show_side_plot)
            if p.type in ["heatmap"]:
                plots_n_cols = 2

        if fig_width and fig_height:
            plt.figure(figsize=[fig_width, fig_height * plots_n_rows], dpi=dpi)
        else:
            plt.figure(dpi=dpi)

        logger.debug(f"plots n_rows={plots_n_rows}; n_cols = {plots_n_cols}")
        logger.info("init graph_coords")
        exon_scale = kwargs.get("exon_scale", 1)
        intron_scale = kwargs.get("intron_scale", .5)
        if plots_n_cols > 1 and intron_scale != 1:
            logger.warning(f"heatmap require intron_scale = 1")
            intron_scale = 1

        self.graph_coords = init_graph_coords(
            self.region, self.exons,
            exon_scale=exon_scale,
            intron_scale=intron_scale
        )

        if plots_n_cols > 1:
            gs = gridspec.GridSpec(plots_n_rows, plots_n_cols, width_ratios=(.99, .01), wspace=0.01, hspace=.15)
        else:
            gs = gridspec.GridSpec(plots_n_rows, plots_n_cols, wspace=.7, hspace=.15)

        curr_idx = 0
        for idx, p in enumerate(self.plots):
            ax_var = plt.subplot(gs[idx, 0])
            if p.type == "density":
                plot_density(
                    ax=ax_var,
                    obj=p.obj[0],
                    graph_coords=self.graph_coords,
                    *args, **kwargs
                )
            elif p.type == "heatmap":
                plot_heatmap(
                    ax=ax_var,
                    cbar_ax=plt.subplot(gs[idx, 1]),
                    data=p.data,
                    y_label=p.group,
                    graph_coords=self.graph_coords,
                    *args, **kwargs
                )
            elif p.type == "line":
                plot_line(
                    ax=ax_var,
                    data=p.data,
                    y_label=p.group,
                    graph_coords=self.graph_coords,
                    *args, **kwargs
                )
            else:
                raise ValueError(f"unknown plot type {p.type}")

            # adjust indicator lines and focus background
            set_indicator_lines(ax=ax_var, sites=self.sites, graph_coords=self.graph_coords)
            set_focus(ax=ax_var, focus=self.focus, graph_coords=self.graph_coords)

            if idx == len(self.plots) - 1:
                set_x_ticks(
                    ax=ax_var, region=self.region,
                    graph_coords=self.graph_coords,
                    sequence=self.sequence.data if self.sequence else None,
                    *args, **kwargs
                )
            curr_idx += 1

        if self.reference is not None:
            ax_var = plt.subplot(gs[curr_idx:curr_idx+self.reference.len(scale=reference_scale), 0])
            plot_reference(ax=ax_var, obj=self.reference, graph_coords=self.graph_coords, *args, **kwargs)

            # adjust indicator lines and focus background
            set_indicator_lines(ax=ax_var, sites=self.sites, graph_coords=self.graph_coords)
            set_focus(ax=ax_var, focus=self.focus, graph_coords=self.graph_coords)
            curr_idx += self.reference.len(scale=reference_scale)

        if self.stroke:
            ax_var = plt.subplot(gs[curr_idx:curr_idx + self.reference.len(scale=reference_scale), 0])
            plot_stroke(ax=ax_var, data=self.stroke, graph_coords=self.graph_coords, *args, **kwargs)

        if output:
            plt.savefig(
                output,
                transparent=True,
                bbox_inches='tight'
            )
        else:
            plt.show()

        return self


if __name__ == '__main__':
    def test_plot():
        from conf.logger import init_logger
        init_logger("INFO")
        Plot().set_reference(
            "../example/example.sorted.gtf.gz"
        ).set_region(
            "chr1", 1270656, 1284730, "+"
        ).add_density(
            path="../example/bams/1.bam",
            category="bam"
        ).add_density(
            path="../example/bws/2.bw",
            category="bw"
        ).add_line(
            path="../example/bams/1.bam",
            category="bam",
            group="1"
        ).add_line(
            path="../example/bams/2.bam",
            category="bam",
            group="1"
        ).add_heatmap(
            path="../example/bams/1.bam",
            category="bam",
            group="1"
        ).add_heatmap(
            path="../example/bams/2.bam",
            category="bam",
            group="1"
        ).add_sites(
            1270656 + 1000
        ).add_sites(
            1270656 + 1000
        ).add_sites(
            1270656 + 2000
        ).add_focus(
            f"{1270656 + 2000}-{1270656 + 3000}"
        ).add_focus(
            f"{1270656 + 5000}-{1270656 + 7000}"
        ).add_stroke(
            f"{1270656 + 5000}-{1270656 + 7000}:{1270656 + 7200}-{1270656 + 8000}@blue"
        ).add_stroke(
            start=1270656 + 7500,
            end=1270656 + 8200,
            color="green",
            label="test"
        ).plot("test_plot.png")

    test_plot()
    pass

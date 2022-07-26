#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""
import math
import os.path

from copy import deepcopy
from typing import List, Optional, Set, Union, Dict

import matplotlib.pyplot as plt
from matplotlib import gridspec

from conf.logger import logger
from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.base.ReadDepth import ReadDepth
from sashimi.file.HiCMatrixTrack import HiCTrack
from sashimi.file.ReadSegments import ReadSegment
from sashimi.base.Stroke import Stroke
from sashimi.file.Bam import Bam
from sashimi.file.Bigwig import Bigwig
from sashimi.file.Depth import Depth
from sashimi.file.Fasta import Fasta
from sashimi.file.File import File
from sashimi.file.Reference import Reference
from sashimi.plot_func import plot_line, plot_density, plot_reference, plot_heatmap, init_graph_coords, set_x_ticks, \
    set_indicator_lines, set_focus, plot_stroke, plot_igv_like, plot_side_plot, plot_hic
from sashimi.file.Junction import load_custom_junction


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

    def __hash__(self) -> int:
        if self.type in ["heatmap", "line"]:
            return hash(self.group)
        return hash((tuple(self.obj), self.group, self.type, tuple(self.category)))

    def __add__(self, other):
        if self.category == ["heatmap", "line"]:
            self.obj[0] += other.obj
        else:
            self.obj += other.obj
        return self

    @property
    def data(self) -> Dict[str, Union[ReadDepth, ReadSegment]]:
        data = {}
        for obj in self.obj:
            if isinstance(obj.data, dict):
                data.update(data)
            elif isinstance(obj, ReadSegment):
                data[obj.label] = obj
            else:
                data[obj.label] = obj.data
        return data

    def len(self, scale: Union[int, float] = .25) -> int:
        n = 0
        if not self.category:
            pass
        elif self.type == "side-plot" and self.category[0] == "bam":
            n += 2
        elif self.type == "igv":
            n += self.obj[0].len(scale / 8)
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

    def load(self, region: GenomicLoci, *args, **kwargs):
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
        self.params = {}
        self.junctions = {}

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

    def set_reference(self, gtf: str,
                      add_domain: bool = False,
                      local_domain: Optional[str] = False,
                      interval: Optional[str] = None,
                      interval_label: Optional[str] = None,
                      transcripts: Optional[List[str]] = None,
                      remove_empty_transcripts: bool = False,
                      color: Optional[str] = None,

                      # transcripts related parameters
                      font_size: int = 5,
                      show_gene: bool = False,
                      show_id: bool = False,
                      reverse_minus: bool = False,
                      exon_width: float = .3,
                      show_exon_id: bool = False,
                      theme: str = "blank"
                      ):
        u"""
        add transcripts to this region
        :param gtf:
        :param add_domain:
        :param local_domain:
        :param interval:
        :param interval_label:
        :param font_size: the size of transcript id, name
        :param transcripts: the list of name or ids of transcripts to draw
        :param remove_empty_transcripts: whether to remove transcripts without any exons
        :param color: the color of exons
        :param show_gene: whether to show gene name/id
        :param show_id: show gene id or gene name
        :param reverse_minus: whether to remove strand of transcripts
        :param theme: the theme of transcript
        :param exon_width: the height of exons
        :param show_exon_id: whether to show exon id
        :return:
        """
        logger.info(f"set reference file to {gtf}")
        self.reference = Reference.create(
            gtf,
            add_domain=add_domain,
            add_local_domain=local_domain
        )

        if interval and interval_label:
            self.reference.add_interval(interval, interval_label)

        self.params["reference"] = {
            "transcripts": transcripts,
            "remove_empty_transcripts": remove_empty_transcripts,
            "color": color,
            "font_size": font_size,
            "show_gene": show_gene,
            "show_id": show_id,
            "reverse_minus": reverse_minus,
            "exon_width": exon_width,
            "show_exon_id": show_exon_id,
            "theme": theme
        }

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
                            barcode_tag: str = "BC",
                            umi_tag: str = "UB",
                            library: str = "fru",
                            features: Optional[dict] = None,
                            deletion_ignore: Optional[int] = True,
                            del_ratio_ignore: float = .5,
                            exon_focus: Optional[str] = None,
                            # for hic plot
                            trans: Optional[str] = None,
                            depth: Optional[int] = 30000
                            ):
        if category == "bam":
            obj = Bam.create(
                path,
                label=label,
                title=title,
                barcodes=barcodes,
                barcode_tag=barcode_tag,
                umi_tag=umi_tag,
                library=library
            )
        elif category == "igv":
            obj = ReadSegment.create(
                path=path,
                label=label,
                library=library,
                features=features,
                deletion_ignore=deletion_ignore,
                del_ratio_ignore=del_ratio_ignore,
                exon_focus=exon_focus
            )
        elif category == "hic":
            obj = HiCTrack.create(
                path=path,
                label=label,
                trans=trans,
                depth=depth
            )

        elif category == "bigwig" or category == "bw":
            category = "bw"
            obj = Bigwig.create(path, label=label, title=title)
        elif category == "depth":
            obj = Depth.create(path, label=label, title=title)
        else:
            raise ValueError(f"the category should be one of [bam, bigwig, bw, depth], instead of {category}")
        return obj, category

    def add_customized_junctions(self, path: str):
        if path and os.path.exists(path) and os.path.isfile(path):
            self.junctions = load_custom_junction(path)

    def add_density(self,
                    path: str,
                    category: str,

                    # file loading parameters
                    label: Union[str, List[str]] = "",
                    title: str = "",
                    barcodes: Optional[Set[str]] = None,
                    barcode_tag: str = "BC",
                    umi_tag: str = "UB",
                    library: str = "fru",

                    # plotting parameters
                    color="blue",
                    font_size: int = 8,
                    show_junction_number: bool = True,
                    junction_number_font_size: int = 5,
                    n_y_ticks: int = 4,
                    distance_between_label_axis: float = .1,
                    show_y_label: bool = True,
                    y_label: str = "",
                    theme: str = "ticks_blank",

                    # side plot parameters
                    show_side_plot: bool = False,
                    strand_choice: Optional[str] = None,
                    ):
        u"""
        add density object to plot
        :param path: the path to input file
        :param category: the input file type
        :param show_side_plot: draw the density distribution of reads from different strand
        :param label: the label of input file
        :param title: the title of input file
        :param barcodes: list of required barcodes
        :param barcode_tag: cell barcode tag
        :param umi_tag: umi barcode tag
        :param library: should be one of [frf: "fr-firststrand", frs:"fr-secondstrand", fru:"fr-unstrand"], default: fru
        :param font_size: the font size for ticks, y-axis label and title
        :param show_junction_number: whether to show the number of junctions
        :param distance_between_label_axis: distance between y-axis label and y-axis ticks
        :param n_y_ticks: number of y ticks
        :param junction_number_font_size:
        :param color: color for this density plot
        :param show_y_label: whether to show y-axis label
        :param y_label: the text of y-axis title
        :param theme: the theme name
        :param strand_choice: the strand to draw on side plot
        :return:
        """
        obj, category = self.__init_input_file__(
            path=path,
            category=category,
            label=label,
            title=title,
            barcodes=barcodes,
            barcode_tag=barcode_tag,
            umi_tag=umi_tag,
            library=library
        )

        type_ = "density"
        if show_side_plot and category == "bam":
            type_ = "side-plot"
        elif show_side_plot:
            logger.warning("show_side_plot only works with bam files")

        info = PlotInfo(obj=obj, type_=type_, category=category)
        self.plots.append(info)
        self.params[info] = {
            "show_junction_number": show_junction_number,
            "junction_number_font_size": junction_number_font_size,
            "color": color,
            "font_size": font_size,
            "n_y_ticks": n_y_ticks,
            "distance_between_label_axis": distance_between_label_axis,
            "show_y_label": show_y_label,
            "y_label": y_label,
            "theme": theme,
            "strand_choice": strand_choice
        }
        return self

    def add_heatmap(self,
                    path: str,
                    group: str,
                    category: str = "bam",

                    # file loading parameters
                    label: Union[str, List[str]] = "",
                    title: str = "",
                    barcodes: Optional[Set[str]] = None,
                    barcode_tag: str = "BC",
                    umi_tag: str = "UB",
                    library: str = "fru",

                    # plotting parameters
                    color="viridis",
                    font_size: int = 8,
                    distance_between_label_axis: float = .1,
                    show_y_label: bool = True,
                    theme: str = "ticks",
                    do_scale: bool = False,
                    clustering: bool = False,
                    clustering_method: str = "ward",
                    distance_metric: str = "euclidean",
                    ):
        u"""
        add multiple objects for a group of heatmap
        :param path: path to input files
        :param group: the heatmap group
        :param category: file category corresponding to input file
        :param label: the label of input file
        :param title: the title of input file
        :param barcodes: list of required barcodes
        :param barcode_tag: cell barcode tag
        :param umi_tag: umi barcode tag
        :param library: fru: fr-unstrand
        :param color: color for this density plot
        :param show_y_label: whether to show y-axis label
        :param theme: the theme name
        :param font_size:
        :param distance_between_label_axis:
        :param do_scale: whether to scale the matrix
        :param clustering: whether reorder matrix by clustering
        :param clustering_method: same as  scipy.cluster.hierarchy.linkage
        :param distance_metric: same as scipy.spatial.distance.pdist
        :param color: used for seaborn.heatmap, see: https://matplotlib.org/3.5.1/tutorials/colors/colormaps.html
                    'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                    'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
                    'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'
        :return:
        """
        obj, category = self.__init_input_file__(
            path=path,
            category=category,
            label=label,
            title=title,
            barcodes=barcodes,
            barcode_tag=barcode_tag,
            umi_tag=umi_tag,
            library=library
        )

        exists = False
        for p in self.plots:
            if p.group == group and p.type == "heatmap":
                p.add(obj=obj, category=category, type_="heatmap")
                exists = True
                break

        if not exists:
            info = PlotInfo(obj=obj, category=category, type_="heatmap", group=group)
            self.plots.append(info)
            self.params[info] = {
                "color": color,
                "font_size": font_size,
                "distance_between_label_axis": distance_between_label_axis,
                "show_y_label": show_y_label,
                "theme": theme,
                "do_scale": do_scale,
                "clustering": clustering,
                "clustering_method": clustering_method,
                "distance_metric": distance_metric
            }

        return self

    def add_line(self,
                 path: str,
                 group: str,
                 category: str = "bam",

                 # file loading parameters
                 label: Union[str, List[str]] = "",
                 title: str = "",
                 barcodes: Optional[Set[str]] = None,
                 barcode_tag: str = "BC",
                 umi_tag: str = "UB",
                 library: str = "fru",

                 # plotting parameters
                 color="blue",
                 font_size: int = 8,
                 distance_between_label_axis: float = .1,
                 show_y_label: bool = True,
                 line_attrs: Optional[Dict] = None,
                 theme: str = "ticks_blank",
                 n_y_ticks: int = 4,
                 show_legend: bool = False,
                 legend_position: str = "upper right",
                 legend_ncol: int = 0
                 ):
        u"""
        add multiple objects for a group of heatmap
        :param path: path to input files
        :param group: the heatmap group
        :param category: file category corresponding to input file
        :param label: the label of input file
        :param title: the title of input file
        :param barcodes: list of required barcodes
        :param barcode_tag: cell barcode tag
        :param umi_tag: umi barcode tag
        :param library: fru: fr-unstrand
        :param distance_between_label_axis: distance between y-axis label and y-axis ticks
        :param n_y_ticks: number of y ticks
        :param color: color for this density plot
        :param show_y_label: whether to show y-axis label
        :param theme: the theme name
        :param font_size:
        :param line_attrs: the additional attributes to control the line, usd by matpltolib.axes.Axes.plot
        :param show_legend: whether to show legend
        :param legend_position:
        :param legend_ncol:
        :return:
        """
        obj, category = self.__init_input_file__(
            path=path,
            category=category,
            label=label,
            title=title,
            barcodes=barcodes,
            barcode_tag=barcode_tag,
            umi_tag=umi_tag,
            library=library
        )

        if not line_attrs:
            line_attrs = {}
        line_attrs["color"] = color

        exists = False
        for p in self.plots:
            if p.group == group and p.type == "line":
                p.add(obj=obj, category=category, type_="line")
                self.params[p]["line_attrs"][obj.label] = line_attrs
                exists = True
                break

        if not exists:
            info = PlotInfo(obj=obj, category=category, type_="line", group=group)
            self.plots.append(info)
            self.params[info] = {
                "line_attrs": {obj.label: line_attrs},
                "show_legend": show_legend,
                "font_size": font_size,
                "n_y_ticks": n_y_ticks,
                "distance_between_label_axis": distance_between_label_axis,
                "show_y_label": show_y_label,
                "theme": theme,
                "legend_position": legend_position,
                "legend_ncol": legend_ncol
            }

        return self

    def add_hic(
            self,
            path: str,
            category: str = "hic",
            label: str = "",
            color: str = "RdYlBu_r",
            trans: Optional[str] = None,
            show_legend: bool = True,
            depth: int = 30000,
            font_size: int = 8,
            n_y_ticks: int = 4,
            distance_between_label_axis: float = .1,
            show_y_label: bool = True,
            theme: str = "ticks"
    ):
        obj, category = self.__init_input_file__(
            path=path,
            category=category,
            label=label,
            depth=depth,
            trans=trans
        )

        info = PlotInfo(obj=obj, category=category, type_="hic")
        self.plots.append(info)
        self.params[info] = {
            "show_legend": show_legend,
            "color": color,
            "y_label": label,
            "font_size": font_size,
            "n_y_ticks": n_y_ticks,
            "distance_between_label_axis": distance_between_label_axis,
            "show_y_label": show_y_label,
            "theme": theme
        }

        return self

    def add_igv(
            self,
            path: str,
            category: str = "igv",
            label: str = "",
            exon_focus: Optional[str] = None,

            # file loading parameters
            library: str = "fru",
            features: Optional[dict] = None,
            deletion_ignore: Optional[int] = True,
            del_ratio_ignore: float = .5,

            # plotting parameters
            exon_color: Optional[str] = None,
            intron_color: Optional[str] = None,
            feature_color: Optional[str] = None,
            exon_width: float = .3,
            font_size: int = 8,
            n_y_ticks: int = 1,
            distance_between_label_axis: float = .1,
            show_y_label: bool = True,
            theme: str = "ticks_blank"
    ):
        u"""
        Add igv-like plot into track
        :param path: path to input files
        :param category: file category for the input file
        :param library: fru: fr-unstrand
        :param features:
        :param exon_focus:
        :param deletion_ignore:
        :param del_ratio_ignore:
        :param label:
        :param exon_color:
        :param intron_color:
        :param feature_color:
        :param exon_width:
        :param font_size:
        :param n_y_ticks:
        :param distance_between_label_axis:
        :param show_y_label:
        :param theme:
        :return:
        """
        obj, category = self.__init_input_file__(
            path=path,
            category=category,
            library=library,
            features=features,
            deletion_ignore=deletion_ignore,
            del_ratio_ignore=del_ratio_ignore,
            exon_focus=exon_focus
        )

        info = PlotInfo(obj=obj, category=category, type_="igv")
        self.plots.append(info)
        self.params[info] = {
            "y_label": label,
            "exon_color": exon_color,
            "intron_color": intron_color,
            "feature_color": feature_color,
            "exon_width": exon_width,
            "font_size": font_size,
            "n_y_ticks": n_y_ticks,
            "distance_between_label_axis": distance_between_label_axis,
            "show_y_label": show_y_label,
            "theme": theme
        }

        return self

    def merge_by_cell(self):

        plots = {}
        for p in self.plots:
            assert isinstance(p, PlotInfo)

            if p.category in ["density", "line"]:
                label = p.obj[0].label
                if label in plots.keys():
                    plots[label] += p
                else:
                    plots[label] = p
            else:
                plots[p.obj[0].label] = p

        self.plots = list(plots.values())

        return self

    def __len__(self) -> int:
        return self.end - self.start + 1

    def copy(self):
        return deepcopy(self)

    def plot(self,
             output: Optional[str] = None,
             reference_scale: Union[int, float] = .25,
             stroke_scale: Union[int, float] = .25,
             dpi: int = 300,
             fig_width: Union[int, float] = 0,
             fig_height: Union[int, float] = 0,
             raster: bool = False,
             *args, **kwargs):
        u"""
        draw image
        :param output: if output is empty then show this image by plt.showfig
        :param reference_scale: to adjust the size of reference plot
        :param stroke_scale: to adjust the size of stroke plot
        :param dpi: the dpi of saved plot
        :param fig_width: the width of figure, if width == 0, the let matplotlib decide the size of image
        :param fig_height: the height of figure, if height == 0, the let matplotlib decide the size of image
        :param raster: plot rasterizer side plot
        """
        assert self.region is not None, f"please set the plotting region first."

        plots_n_rows = 1
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
                p.load(region=self.region, junctions=self.junctions.get(p.obj[0].label, {}), *args, **kwargs)
            except Exception as err:
                logger.warning(f"failed to load data from {p}")
                raise err

            plots_n_rows += p.len(reference_scale)
            if p.type in ["heatmap", "hic"]:
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

        max_used_y_val = None
        if kwargs.get("same_y"):
            for p in self.plots:
                if p.type in ["density", "side-plot", "line"]:
                    for obj in p.obj:
                        y = max(obj.data.wiggle)

                        if obj.log_trans == "2":
                            y = math.log2(y)
                        elif obj.log_trans == "10":
                            y = math.log10(y)
                        elif obj.log_trans == "1p":
                            y = math.log1p(y)
                        elif obj.log_trans == "e":
                            y = math.log(y)

                        max_used_y_val = y if max_used_y_val is None else max(y, max_used_y_val)

        curr_idx = 0
        for p in self.plots:
            if p.type == "igv":
                ax_var = plt.subplot(gs[curr_idx: curr_idx + p.len(reference_scale), 0])
            else:
                ax_var = plt.subplot(gs[curr_idx, 0])
            if p.type == "density":
                plot_density(
                    ax=ax_var,
                    obj=p.obj[0],
                    graph_coords=self.graph_coords,
                    max_used_y_val=max_used_y_val,
                    **self.params[p]
                )
            elif p.type == "hic":
                plot_hic(
                    ax=ax_var,
                    cbar_ax=plt.subplot(gs[curr_idx, 1]),
                    obj=p.obj,
                    **self.params[p]
                )
            elif p.type == "side-plot":
                plot_density(
                    ax=ax_var,
                    obj=p.obj[0],
                    graph_coords=self.graph_coords,
                    max_used_y_val=max_used_y_val,
                    **self.params[p]
                )

                side_ax = plt.subplot(gs[curr_idx + 1, 0])

                plot_side_plot(
                    side_ax, p.obj[0],
                    graph_coords=self.graph_coords,
                    raster=raster,
                    **self.params[p]
                )
                curr_idx += 1
            elif p.type == "heatmap":
                plot_heatmap(
                    ax=ax_var,
                    cbar_ax=plt.subplot(gs[curr_idx, 1]),
                    data=p.data,
                    y_label=p.group,
                    graph_coords=self.graph_coords,
                    raster=raster,
                    **self.params[p]
                )
            elif p.type == "line":
                plot_line(
                    ax=ax_var,
                    data=p.data,
                    y_label=p.group,
                    max_used_y_val=max_used_y_val,
                    graph_coords=self.graph_coords,
                    **self.params[p]
                )
            elif p.type == "igv":
                plot_igv_like(
                    ax=ax_var,
                    obj=p.data,
                    graph_coords=self.graph_coords,
                    raster=raster,
                    **self.params[p]
                )
            else:
                raise ValueError(f"unknown plot type {p.type}")

            # adjust indicator lines and focus background
            set_indicator_lines(ax=ax_var, sites=self.sites, graph_coords=self.graph_coords)
            set_focus(ax=ax_var, focus=self.focus, graph_coords=self.graph_coords)

            if p.type != "igv":
                curr_idx += 1
            else:
                curr_idx += p.len(reference_scale)

        # draw x label
        set_x_ticks(
            ax=plt.subplot(gs[curr_idx, 0]), region=self.region,
            graph_coords=self.graph_coords,
            sequence=self.sequence.data if self.sequence else None,
            *args, **kwargs
        )
        curr_idx += 1

        if self.reference is not None:
            ax_var = plt.subplot(gs[curr_idx:curr_idx + self.reference.len(scale=reference_scale), 0])
            plot_reference(ax=ax_var, obj=self.reference,
                           graph_coords=self.graph_coords,
                           plot_domain=self.reference.add_domain,
                           **self.params["reference"])

            # adjust indicator lines and focus background
            set_indicator_lines(ax=ax_var, sites=self.sites, graph_coords=self.graph_coords)
            set_focus(ax=ax_var, focus=self.focus, graph_coords=self.graph_coords)
            curr_idx += self.reference.len(scale=reference_scale)

        if self.stroke:
            ax_var = plt.subplot(gs[curr_idx:curr_idx + self.reference.len(scale=reference_scale), 0])
            plot_stroke(ax=ax_var, data=self.stroke, graph_coords=self.graph_coords, *args, **kwargs)

        if output:
            plt.savefig(output, transparent=True, bbox_inches='tight')
        else:
            plt.show()

        plt.close()
        return self


if __name__ == '__main__':
    def test_plot():
        from conf.logger import init_logger
        init_logger("INFO")
        Plot().set_reference(
            "../example/example.sorted.gtf.gz",
            add_domain=True,
            interval="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.bed.gz",
            interval_label="polyA",
            show_gene=True,
            color="pink",
        ).add_interval(
            interval="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.simple.bed.gz",
            interval_label="polyAS"
        ).set_region(
            "chr1", 1270656, 1284730, "+"
        ).add_density(
            path="../example/bams/1.bam",
            category="bam",
            color="blue",
            show_side_plot=True,
        ).add_density(
            path="../example/bws/2.bw",
            category="bw",
            color="green"
        ).add_line(
            path="../example/bams/1.bam",
            category="bam",
            group="1",
            line_attrs={"lw": 3}
        ).add_line(
            path="../example/bams/2.bam",
            category="bam",
            group="2",
            color="red",
            line_attrs={"linestyle": "dashed"}
        ).add_heatmap(
            path="../example/bams/1.bam",
            category="bam",
            group="1",
        ).add_heatmap(
            path="../example/bams/2.bam",
            category="bam",
            group="1"
        ).add_hic(
            path="../example/Li_et_al_2015.h5",
            category="hic",
            trans="log2",
            depth=30000
        ).add_igv(
            path="../example/bams/3.bam",
            features={
                "m6a": "ma",
                "real_strand": "rs",
                "polya": "pa"
            },
            category="igv",
            label="igv"
        ).add_igv(
            path="../example/SRX9697989.corrected_reads.bed.gz",
            category="igv",
            label="bed12"
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
        ).plot("test_plot.png", fig_width=6, fig_height=2, raster=True)


    test_plot()
    pass

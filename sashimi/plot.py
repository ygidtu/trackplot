#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created by ygidtu@gmail.com at 2019.12.06
"""
import copy
import faulthandler
import io
import logging
import os
from multiprocessing import Pool
from typing import Set

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import FigureCanvasPdf

from sashimi.file.ATAC import ATAC
from sashimi.file.Bam import Bam
from sashimi.file.BedGraph import Bedgraph
from sashimi.file.Bigwig import Bigwig
from sashimi.file.Depth import Depth
from sashimi.file.Fasta import Fasta
from sashimi.file.Junction import load_custom_junction
from sashimi.file.Motif import Motif
from sashimi.plot_func import *

logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

faulthandler.enable()


def __load__(args):
    p, args, kwargs = args
    p.load(*args, **kwargs)
    return p


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
        return hash((tuple([hash(x) for x in self.obj]), self.group, self.type, tuple(self.category)))

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __ne__(self, other):
        return hash(self) != hash(other)

    def __add__(self, other):
        if self.category == ["heatmap", "line"]:
            self.obj[0] += other.obj
        else:
            self.obj += other.obj
        return self

    @property
    def is_single_cell(self):
        try:
            return self.obj[0].is_single_cell
        except AttributeError:
            return False

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
        elif self.type == "site-plot" and self.category[0] == "bam":
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

    def load(self, region: GenomicLoci, n_jobs: int = 0, *args, **kwargs):
        if n_jobs <= 1:
            for obj in self.obj:
                obj.load(region=region, *args, **kwargs)
        else:
            with Pool(n_jobs) as p:
                kwargs_ = deepcopy(kwargs)
                kwargs_["region"] = region
                self.obj = p.map(__load__, [[o, args, kwargs_] for o in self.obj])

        if self.type == "density":
            if len(self.obj) > 1:
                data = self.obj[0].data
                for obj in self.obj[1:]:
                    data += obj.data

                for obj in self.obj:
                    obj.data = data

        return self


class Plot(object):
    u"""
    this class is the main framework of sashimi
    """

    def __init__(self):
        u"""
        init this class
        """
        self.__n_objs__ = 0
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
        self.link = []

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

    def set_region(self, chromosome: str = "",
                   start: int = 0, end: int = 0,
                   strand: str = "+",
                   region: Optional[GenomicLoci] = None):
        u"""
        change the plot region
        """

        if region:
            self.region = region
        elif chromosome and start and end:
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

    def add_links(
            self,
            link: Optional[str] = None,
            start: int = 0,
            end: int = 0,
            label: str = "",
            color: str = "blue"
    ):
        u"""
        add stroke to plot
        :param link: string in 100-200@red;300-400 format
        :param start: start site of stroke
        :param end: end site of stroke
        :param label: label of stroke
        :param color: color of stroke
        :return:
        """
        assert self.region is not None, f"please set plot region first."
        if link:
            self.link.append(Stroke.create(link, self.region, default_color=color))

        if 0 < start < end:
            self.link.append([Stroke(start - self.start, end - self.end, color, label)])
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
                      choose_primary: bool = False,
                      color: Optional[str] = "black",

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
        :param choose_primary: whether to choose primary transcript for each gene
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
            "choose_primary": choose_primary,
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

    def __init_input_file__(self, path: str,
                            category: str = "bam",
                            label: Union[str, List[str]] = "",
                            title: str = "",
                            barcode: str = "",
                            barcode_groups: Dict[str, Set[str]] = None,
                            barcode_tag: str = "BC",
                            umi_tag: str = "UB",
                            library: str = "fru",
                            features: Optional[dict] = None,
                            deletion_ignore: Optional[int] = True,
                            del_ratio_ignore: float = .5,
                            exon_focus: Optional[str] = None,
                            density_by_strand: bool = False,

                            # for hic plot
                            trans: Optional[str] = None,
                            depth: Optional[int] = 30000,

                            # for ATAC
                            size_factor=None,

                            log_trans: Optional[str] = None,
                            ):
        self.__n_objs__ += 1
        path = os.path.expanduser(path)
        if category == "bam":
            obj = Bam.create(
                path,
                label=label,
                title=title,
                barcodes=barcode_groups.get(barcode) if barcode_groups else None,
                barcode_tag=barcode_tag,
                umi_tag=umi_tag,
                library=library,
                density_by_strand=density_by_strand
            )
        elif category == "atac":
            obj = ATAC.create(
                path,
                label=label,
                title=title,
                barcode=barcode,
                barcode_groups=barcode_groups,
                size_factors=size_factor
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
        elif category == "bedgraph" or category == "bg":
            category = "bg"
            obj = Bedgraph.create(path=path, label=label, title=title)
        elif category == "depth":
            obj = Depth.create(path, label=label, title=title)
        else:
            raise ValueError(
                f"the category should be one of [bam, bigwig, bw, depth, bedgraph, bg], instead of {category}")

        obj.log_trans = log_trans
        return obj, category

    def add_customized_junctions(self, path: str):
        if path and os.path.exists(path) and os.path.isfile(path):
            self.junctions = load_custom_junction(path)

    def add_density(self,
                    path: str,
                    category: str = "bam",
                    size_factor=None,

                    # file loading parameters
                    label: Union[str, List[str]] = "",
                    title: str = "",
                    barcode: str = "",
                    barcode_groups: Dict[str, Set[str]] = None,
                    barcode_tag: str = "BC",
                    umi_tag: str = "UB",
                    library: str = "fru",
                    density_by_strand: bool = False,

                    # plotting parameters
                    color="blue",
                    font_size: int = 8,
                    show_junction_number: bool = True,
                    junction_number_font_size: int = 5,
                    n_y_ticks: int = 4,
                    show_y_label: bool = True,
                    y_label: str = "",
                    theme: str = "ticks_blank",
                    log_trans: Optional[str] = None,

                    # site plot parameters
                    show_site_plot: bool = False,
                    strand_choice: Optional[str] = None,

                    only_customized_junction: bool = False
                    ):
        u"""
        add density object to plot
        :param path: the path to input file
        :param category: the input file type
        :param show_site_plot: draw the density distribution of reads from different strand
        :param label: the label of input file
        :param title: the title of input file
        :param size_factor
        :param barcode: key of barcode barcode_groups
        :param barcode_groups:
        :param barcode_tag: cell barcode tag
        :param umi_tag: umi barcode tag
        :param density_by_strand: whether to draw density plot in strand-specific manner.
        :param library: should be one of [frf: "fr-firststrand", frs:"fr-secondstrand", fru:"fr-unstrand"], default: fru
        :param font_size: the font size for ticks, y-axis label and title
        :param show_junction_number: whether to show the number of junctions
        :param n_y_ticks: number of y ticks
        :param junction_number_font_size:
        :param color: color for this density plot
        :param show_y_label: whether to show y-axis label
        :param y_label: the text of y-axis title
        :param theme: the theme name
        :param strand_choice: the strand to draw on site plot
        :param only_customized_junction: only draw customized junctions
        :return:
        """
        obj, category = self.__init_input_file__(
            path=path,
            category=category,
            label=label,
            title=title,
            barcode=barcode,
            barcode_groups=barcode_groups,
            barcode_tag=barcode_tag,
            umi_tag=umi_tag,
            library=library,
            size_factor=size_factor,
            density_by_strand=density_by_strand,
            log_trans=log_trans
        )

        type_ = "density"
        if show_site_plot and category == "bam":
            type_ = "site-plot"
        elif show_site_plot:
            logger.warning("show_site_plot only works with bam files")

        info = PlotInfo(obj=obj, type_=type_, category=category)

        new_obj = True
        for p in self.plots:
            if p.category == info.category:
                for obj_ in p.obj:
                    if obj_.label == label and label:
                        param = self.params.pop(p)
                        p.obj.append(obj)
                        new_obj = False
                        self.params[p] = param
                        break
            if not new_obj:
                break

        if new_obj:
            self.plots.append(info)
            self.params[info] = {
                "show_junction_number": show_junction_number,
                "junction_number_font_size": junction_number_font_size,
                "color": color,
                "font_size": font_size,
                "n_y_ticks": n_y_ticks,
                "show_y_label": show_y_label,
                "y_label": y_label,
                "theme": theme,
                "strand_choice": strand_choice,
                "density_by_strand": density_by_strand,
                "only_customized_junction": only_customized_junction
            }
        return self

    def add_heatmap(self,
                    path: str,
                    group: str = "",
                    category: str = "bam",
                    size_factor=None,

                    # file loading parameters
                    label: Union[str, List[str]] = "",
                    title: str = "",
                    barcode: str = "",
                    barcode_groups: Dict[str, Set[str]] = None,
                    barcode_tag: str = "BC",
                    umi_tag: str = "UB",
                    library: str = "fru",

                    # plotting parameters
                    color="viridis",
                    font_size: int = 8,
                    show_y_label: bool = True,
                    theme: str = "ticks_blank",
                    do_scale: bool = False,
                    clustering: bool = False,
                    clustering_method: str = "ward",
                    distance_metric: str = "euclidean",
                    show_row_names: bool = False,
                    vmin=None, vmax=None,
                    log_trans: Optional[str] = None,
                    ):
        u"""
        add multiple objects for a group of heatmap
        :param path: path to input files
        :param group: the heatmap group
        :param category: file category corresponding to input file
        :param label: the label of input file
        :param size_factor: only used by atac
        :param title: the title of input file
        :param barcode: key of barcode barcode_groups
        :param barcode_groups:
        :param barcode_tag: cell barcode tag
        :param umi_tag: umi barcode tag
        :param library: fru: fr-unstrand
        :param color: color for this density plot
        :param show_y_label: whether to show y-axis label
        :param theme: the theme name
        :param font_size:
        :param do_scale: whether to scale the matrix
        :param clustering: whether reorder matrix by clustering
        :param clustering_method: same as  scipy.cluster.hierarchy.linkage
        :param distance_metric: same as scipy.spatial.distance.pdist
        :param color: used for seaborn.heatmap, see: https://matplotlib.org/3.5.1/tutorials/colors/colormaps.html
                    'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
                    'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
                    'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'
        :param show_row_names:
        :param vmin: Values to anchor the colormap,
                     otherwise they are inferred from the data and other keyword arguments.
        :param vmax: Values to anchor the colormap,
                     otherwise they are inferred from the data and other keyword arguments.
        :return:
        """
        obj, category = self.__init_input_file__(
            path=path,
            category=category,
            label=label,
            title=title,
            barcode=barcode,
            barcode_groups=barcode_groups,
            barcode_tag=barcode_tag,
            umi_tag=umi_tag,
            library=library,
            size_factor=size_factor,
            log_trans=log_trans
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
                "show_y_label": show_y_label,
                "theme": theme,
                "do_scale": do_scale,
                "clustering": clustering,
                "clustering_method": clustering_method,
                "distance_metric": distance_metric,
                "show_row_names": show_row_names,
                "vmin": vmin, "vmax": vmax
            }

        return self

    def add_line(self,
                 path: str,
                 group: str = "",
                 category: str = "bam",

                 # file loading parameters
                 label: Union[str, List[str]] = "",
                 title: str = "",
                 barcode: str = "",
                 barcode_groups: Dict[str, Set[str]] = None,
                 barcode_tag: str = "BC",
                 umi_tag: str = "UB",
                 library: str = "fru",

                 # plotting parameters
                 color="blue",
                 font_size: int = 8,
                 show_y_label: bool = True,
                 line_attrs: Optional[Dict] = None,
                 theme: str = "ticks_blank",
                 n_y_ticks: int = 4,
                 show_legend: bool = False,
                 legend_position: str = "upper right",
                 legend_ncol: int = 0,
                 log_trans: Optional[str] = None,
                 ):
        u"""
        add multiple objects for a group of heatmap
        :param path: path to input files
        :param group: the heatmap group
        :param category: file category corresponding to input file
        :param label: the label of input file
        :param title: the title of input file
        :param barcode: key of barcode barcode_groups
        :param barcode_groups:
        :param barcode_tag: cell barcode tag
        :param umi_tag: umi barcode tag
        :param library: fru: fr-unstrand
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
            barcode=barcode,
            barcode_groups=barcode_groups,
            barcode_tag=barcode_tag,
            umi_tag=umi_tag,
            library=library,
            log_trans=log_trans
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
            "show_y_label": show_y_label,
            "theme": theme
        }

        return self

    def add_motif(self,
                  path: str,
                  category: str = "motif",
                  motif_region: GenomicLoci = None,
                  width: float = 0.8,
                  theme: str = "blank",
                  **kwargs):
        obj = Motif.create(path, self.region)

        info = PlotInfo(obj=obj, category=category, type_="motif")

        self.plots.append(info)
        self.params[info] = {"width": width, "theme": theme, "region": motif_region if motif_region else self.region}
        return self

    def add_manual(self,
                   data: np.array,
                   image_type: str = "line",
                   label: str = "",
                   group: str = "",
                   color: str = "blue",
                   font_size: int = 8,
                   n_y_ticks: int = 1,
                   show_y_label: bool = True,
                   theme: str = "ticks_blank",
                   **kwargs
                   ):
        obj = File("")
        obj.label = label
        obj.data = ReadDepth(data)
        info = PlotInfo(obj=obj, category="manual", type_=image_type, group=group)

        exists = False
        if group:
            for p in self.plots:
                if p.group == group and p.type == image_type:
                    p.add(obj=obj, category="manual", type_=image_type)
                    self.params[p]["line_attrs"][obj.label] = {"color": color}
                    exists = True
                    break

        if not exists:
            self.plots.append(info)
            self.params[info] = {
                "font_size": font_size,
                "n_y_ticks": n_y_ticks,
                "show_y_label": show_y_label,
                "theme": theme,
                "color": color,
                "line_attrs": {obj.label: {"color": color}}
            }
            self.params[info].update(kwargs)

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
             width: Union[int, float] = 0,
             height: Union[int, float] = 0,
             raster: bool = False,
             return_image: Optional[str] = None,
             sc_height_ratio: Optional[Dict[str, float]] = None,
             distance_between_label_axis: float = .3,
             n_jobs: int = 1,
             *args, **kwargs):
        u"""
        draw image
        :param output: if output is empty then show this image by plt.showfig
        :param reference_scale: to adjust the size of reference plot
        :param stroke_scale: to adjust the size of stroke plot
        :param dpi: the dpi of saved plot
        :param width: the width of figure, if width == 0, the let matplotlib decide the size of image
        :param height: the height of figure, if height == 0, the let matplotlib decide the size of image
        :param raster: plot rasterizer site plot
        :param sc_height_ratio: adjust the relative height of single cell plots
        :param distance_between_label_axis: distance between y-axis label and y-axis ticks
        :param n_jobs: load data in how many processes
        :param return_image: used for interactive ui
        """
        if sc_height_ratio is None:
            sc_height_ratio = {"density": .2, "heatmap": .2}
        assert self.region is not None, f"please set the plotting region first."

        plots_n_rows = 1
        plots_n_cols = 1

        height_ratio = []
        if self.reference is not None:
            logger.info("load reference")
            self.reference.load(self.region, *args, **kwargs)
            plots_n_rows += self.reference.len(scale=reference_scale)

        if self.stroke:
            plots_n_rows += int(max(len(self.stroke) * stroke_scale, 1))

        if self.link:
            plots_n_rows += len(self.link)

        if self.sequence is not None:
            logger.info("load sequence")
            self.sequence.load(self.region)

        logger.info(f"load data of {len(self.plots)} plots")

        cmds = []
        for p in self.plots:
            assert isinstance(p, PlotInfo), f"unrecognized data type: {type(p)}"
            # only enable the multiprocessing while n_jobs > 1
            if self.__n_objs__ / len(self.plots) >= n_jobs > 1:
                p.load(self.region, n_jobs, junctions=self.junctions.get(p.obj[0].label, {}), *args, **kwargs)
            elif n_jobs > 1:
                temp = copy.deepcopy(kwargs)
                temp["region"] = self.region if p.type != "motif" else self.params[p]["region"]
                temp["junctions"] = self.junctions.get(p.obj[0].label, {})
                cmds.append([p, args, temp])

        if len(cmds) > 0:
            with Pool(max(1, min(n_jobs, len(self.plots)))) as p:
                self.plots = p.map(__load__, cmds)

        # count the plots size
        for p in self.plots:
            if n_jobs <= 1:
                p.load(self.region, junctions=self.junctions.get(p.obj[0].label, {}), *args, **kwargs)

            for obj in p.obj:
                obj.transform()

            plots_n_rows += p.len(reference_scale)
            if p.type in ["heatmap", "hic"]:
                plots_n_cols = 2

            if p.is_single_cell:
                height_ratio.append(sc_height_ratio.get(p.type, 1))
            else:
                height_ratio.append(1)

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

        height_ratio += [1 for _ in range(plots_n_rows - len(height_ratio))]

        if width and height:
            fig = plt.figure(figsize=[width, height * sum(height_ratio)], dpi=dpi)
        else:
            fig = plt.figure(dpi=dpi)

        if plots_n_cols > 1:
            gs = gridspec.GridSpec(plots_n_rows, plots_n_cols, height_ratios=height_ratio,
                                   width_ratios=(.99, .01), wspace=0.01, hspace=.15)
        else:
            gs = gridspec.GridSpec(plots_n_rows, plots_n_cols, height_ratios=height_ratio,
                                   wspace=.7, hspace=.15)

        max_used_y_val = None
        if kwargs.get("same_y"):
            for p in self.plots:
                if p.type in ["density", "site-plot", "line"]:
                    for obj in p.obj:
                        y = max(obj.data.wiggle)
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
                    distance_between_label_axis=distance_between_label_axis,
                    raster=raster,
                    **self.params[p]
                )
            elif p.type == "hic":
                plot_hic(
                    ax=ax_var,
                    cbar_ax=plt.subplot(gs[curr_idx, 1]),
                    obj=p.obj,
                    distance_between_label_axis=distance_between_label_axis,
                    raster=raster,
                    **self.params[p]
                )
            elif p.type == "site-plot":
                plot_density(
                    ax=ax_var,
                    obj=p.obj[0],
                    graph_coords=self.graph_coords,
                    max_used_y_val=max_used_y_val,
                    distance_between_label_axis=distance_between_label_axis,
                    raster=raster,
                    **self.params[p]
                )

                site_ax = plt.subplot(gs[curr_idx + 1, 0])

                plot_site_plot(
                    site_ax, p.obj[0],
                    graph_coords=self.graph_coords,
                    raster=raster,
                    distance_between_label_axis=distance_between_label_axis,
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
                    distance_between_label_axis=distance_between_label_axis,
                    **self.params[p]
                )
            elif p.type == "line":
                plot_line(
                    ax=ax_var,
                    data=p.data,
                    y_label=p.group,
                    max_used_y_val=max_used_y_val,
                    graph_coords=self.graph_coords,
                    distance_between_label_axis=distance_between_label_axis,
                    **self.params[p]
                )
            elif p.type == "igv":
                plot_igv_like(
                    ax=ax_var,
                    obj=p.data,
                    graph_coords=self.graph_coords,
                    raster=raster,
                    distance_between_label_axis=distance_between_label_axis,
                    **self.params[p]
                )
            elif p.type == "motif":
                plot_motif(ax=ax_var, obj=p.obj[0], graph_coords=self.graph_coords, **self.params[p])
            else:
                raise ValueError(f"unknown plot type {p.type}")

            # adjust indicator lines and focus background
            set_indicator_lines(ax=ax_var, sites=self.sites, graph_coords=self.graph_coords)
            set_focus(ax=ax_var, focus=self.focus, graph_coords=self.graph_coords)

            if p.type != "igv":
                curr_idx += 1
            else:
                curr_idx += p.len(reference_scale)

        if self.link:
            for link in self.link:
                plot_links(ax=plt.subplot(gs[curr_idx:(curr_idx + 1), 0]),
                           data=link, graph_coords=self.graph_coords)
                curr_idx += 1

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
                           distance_between_label_axis=distance_between_label_axis,
                           **self.params["reference"])

            # adjust indicator lines and focus background
            set_indicator_lines(ax=ax_var, sites=self.sites, graph_coords=self.graph_coords)
            set_focus(ax=ax_var, focus=self.focus, graph_coords=self.graph_coords)
            curr_idx += self.reference.len(scale=reference_scale)

        if self.stroke:
            ax_var = plt.subplot(gs[curr_idx:plots_n_rows, 0])
            plot_stroke(ax=ax_var, data=self.stroke, graph_coords=self.graph_coords, *args, **kwargs)

        if output:
            fig.savefig(output, transparent=True, bbox_inches='tight')
        elif return_image:
            output = io.BytesIO()
            if return_image == "png":
                FigureCanvasAgg(fig).print_png(output)
            elif return_image == "pdf":
                FigureCanvasPdf(fig).print_pdf(output)
            return output
        else:
            plt.show()

        plt.close()
        return self


if __name__ == '__main__':
    pass

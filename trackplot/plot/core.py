#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
Plot — the main orchestration class that manages configuration, data loading,
layout computation, and plot assembly.

Dependencies: all other plot/ submodules, file modules
"""

import copy
import io
import logging
import os
from datetime import datetime
from multiprocessing import Pool
from typing import List, Optional

from loguru import logger
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import FigureCanvasPdf

from trackplot.base.GenomicLoci import GenomicLoci
from trackplot.base.ReadDepth import ReadDepth
from trackplot.base.Stroke import Stroke
from trackplot.file.Annotation import Annotation
from trackplot.file.Bam import Bam
from trackplot.file.BedGraph import Bedgraph
from trackplot.file.Bigwig import Bigwig
from trackplot.file.Depth import Depth
from trackplot.file.Fasta import Fasta
from trackplot.file.File import File
from trackplot.file.HiCMatrixTrack import HiCTrack
from trackplot.file.Junction import load_custom_junction
from trackplot.file.Motif import Motif
from trackplot.file.ReadSegments import ReadSegment
from trackplot.plot.coord import (init_graph_coords, set_focus,
                                  set_indicator_lines, set_x_ticks)
from trackplot.plot.info import PlotInfo, __load__
from trackplot.plot.limits import precompute_y_limits
from trackplot.plot.render import (plot_annotation, plot_density, plot_heatmap,
                                   plot_hic, plot_igv_like, plot_line,
                                   plot_links, plot_motif, plot_site_plot,
                                   plot_stroke)

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


def _resolve_version():
    """Read version from installed package metadata, falling back to pyproject.toml."""
    try:
        from importlib.metadata import version

        return version("trackplot")
    except Exception:
        pass
    try:
        import tomllib

        _here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(_here, "..", "..", "pyproject.toml"), "rb") as f:
            return tomllib.load(f)["project"]["version"]
    except Exception:
        return "0.0.0"


__version__ = _resolve_version()
__author__ = "Yiming Zhang & Ran Zhou"
__email__ = "ygidtu@gmail.com"


def add_object_error(func):
    """Decorator: catch and log errors from plot configuration methods."""

    def inner(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as err:
            logger.debug(f"trackplot will ignore this object, {err}")

    return inner


class Plot(object):
    """
    This class is the main framework of sashimi plot.
    It manages region, annotation, file tracks, layout, and rendering.
    """

    __slots__ = [
        "__n_objs__",
        "region",
        "sites",
        "focus",
        "stroke",
        "events",
        "sequence",
        "annotation",
        "graph_coords",
        "plots",
        "params",
        "junctions",
        "link",
    ]

    def __init__(
        self,
        logfile: Optional[str] = None,
        backend: Optional[str] = None,
        font_family: Optional[str] = None,
    ):
        self.__n_objs__ = 0
        self.region = None
        self.sites = {}
        self.focus = {}
        self.stroke = []
        self.events = None
        self.sequence = None
        self.annotation = None
        self.graph_coords = None
        self.plots = []
        self.params = {}
        self.junctions = {}
        self.link = []

        if logfile:
            logger.add(logfile, level="DEBUG")
        logger.info(f"Create trackplot version: {__version__}")

        if backend:
            try:
                plt.switch_backend(backend)
            except ImportError as err:
                if backend.lower() == "cairo":
                    logger.debug("Cairo backend required cairocffi installed")
                    logger.debug("Switch back to Agg backend")
                else:
                    logger.debug(f"backend error, switch back to Agg: {err}")
                plt.switch_backend("Agg")

        plt.rcParams["pdf.fonttype"] = 42
        if font_family:
            plt.rcParams["font.family"] = font_family

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def chrom(self) -> Optional[str]:
        return self.region.chrom if self.region else None

    @property
    def start(self) -> Optional[int]:
        return self.region.start if self.region else None

    @property
    def end(self) -> Optional[int]:
        return self.region.end if self.region else None

    @property
    def strand(self) -> Optional[str]:
        return self.region.strand if self.region else None

    @property
    def exons(self) -> Optional[List[List[int]]]:
        return self.annotation.exons if self.annotation else None

    # ------------------------------------------------------------------
    # Configuration methods
    # ------------------------------------------------------------------

    @add_object_error
    def set_region(
        self,
        chromosome: str = "",
        start: int = 0,
        end: int = 0,
        strand: str = "+",
        region: Optional[GenomicLoci] = None,
    ):
        if region:
            self.region = region
        elif chromosome and start and end:
            self.region = GenomicLoci(chromosome, start=start, end=end, strand=strand)
        logger.info(f"set region to {self.region}")
        return self

    @add_object_error
    def add_sites(self, sites):
        assert self.region is not None, "please set plot region first."
        if sites is not None:
            logger.info(f"add sites: {sites}")
            if isinstance(sites, str):
                sites = [int(x) for x in sites.split(",")]
                for s in sites:
                    s = s - self.start
                    self.sites[s] = "red" if s in self.sites else "blue"
            elif isinstance(sites, int):
                sites = sites - self.start
                self.sites[sites] = "red" if sites in self.sites else "blue"
        return self

    @add_object_error
    def add_focus(self, focus: Optional[str] = None, start: int = 0, end: int = 0):
        assert self.region is not None, "please set plot region first."
        if focus:
            logger.info(f"add focus: {focus}")
            for site in focus.split(":"):
                site = sorted([int(x) - self.start for x in site.split("-")])
                if site[0] < 0:
                    site[0] = 0
                if site[-1] > len(self):
                    site[-1] = len(self)
                self.focus[site[0]] = max(site[1], self.focus.get(site[0], -1))
        if 0 < start < end:
            logger.info(f"add focus: {start}-{end}")
            self.focus[start] = max(end, self.focus.get(start, -1))
        return self

    @add_object_error
    def add_stroke(
        self,
        stroke: Optional[str] = None,
        start: int = 0,
        end: int = 0,
        label: str = "",
        color: str = "black",
    ):
        assert self.region is not None, "please set plot region first."
        if stroke:
            logger.info(f"add stroke: {stroke}")
            self.stroke += Stroke.create(stroke, self.region)
        if 0 < start < end:
            logger.info(f"add stroke: {start}-{end}")
            self.stroke.append(Stroke(start - self.start, end - self.end, color, label))
        return self

    @add_object_error
    def add_links(
        self,
        link: Optional[str] = None,
        start: int = 0,
        end: int = 0,
        label: str = "",
        color: str = "blue",
    ):
        assert self.region is not None, "please set plot region first."
        if link:
            logger.info(f"add link: {link}")
            self.link.append(Stroke.create(link, self.region, default_color=color))
        if 0 < start < end:
            logger.info(f"add link: {start}-{end}")
            self.link.append([Stroke(start - self.start, end - self.end, color, label)])
        return self

    @add_object_error
    def set_sequence(self, fasta: str):
        logger.info(f"fetch sequence from {fasta}")
        self.sequence = Fasta.create(fasta)
        return self

    @add_object_error
    def set_annotation(self, gtf: str, **kwargs):
        logger.info(f"set annotation file to {gtf}")
        self.annotation = Annotation.create(
            gtf,
            add_domain=kwargs.get("add_domain", False),
            add_local_domain=kwargs.get("local_domain", False),
            domain_exclude=kwargs.get("domain_exclude", False),
            domain_include=kwargs.get("domain_include", False),
        )
        if kwargs.get("interval") and kwargs.get("interval_label"):
            self.add_interval(kwargs["interval"], kwargs["interval_label"])

        self.params["annotation"] = {
            "transcripts": kwargs.get("transcripts"),
            "remove_empty_transcripts": kwargs.get("remove_empty_transcripts", False),
            "choose_primary": kwargs.get("choose_primary", False),
            "color": kwargs.get("color", "black"),
            "font_size": kwargs.get("font_size", 5),
            "show_gene": kwargs.get("show_gene", True),
            "show_id": kwargs.get("show_id", False),
            "exon_width": kwargs.get("exon_width", 0.3),
            "show_exon_id": kwargs.get("show_exon_id", False),
            "theme": kwargs.get("theme", "blank"),
        }
        return self

    @add_object_error
    def add_interval(self, interval: str, interval_label: str):
        assert self.annotation is not None, "please set_annotation first."
        logger.info(f"add interval: {interval} - {interval_label}")
        self.annotation.add_interval(interval, interval_label)
        return self

    # ------------------------------------------------------------------
    # Input file initialization
    # ------------------------------------------------------------------

    def __init_input_file__(
        self,
        path,
        category="bam",
        label="",
        title="",
        barcode="",
        barcode_groups=None,
        barcode_tag="BC",
        umi_tag="UB",
        library="fru",
        features=None,
        deletion_ignore=True,
        del_ratio_ignore=0.5,
        exon_focus=None,
        density_by_strand=False,
        depth=30000,
        tad=None,
        size_factor=None,
        log_trans=None,
        **kwargs,
    ):
        self.__n_objs__ += 1
        path = os.path.expanduser(path)
        logger.info(f"add {category} {label} {path}")

        if category == "bam":
            obj = Bam.create(
                path,
                label=label,
                title=title,
                barcodes=barcode_groups.get(barcode) if barcode_groups else None,
                barcode_tag=barcode_tag,
                umi_tag=umi_tag,
                library=library,
                density_by_strand=density_by_strand,
                size_factor=size_factor,
            )
        elif category == "atac":
            from trackplot.file.ATAC import ATAC

            obj = ATAC.create(
                path,
                label=label,
                title=title,
                barcode=barcode,
                barcode_groups=barcode_groups,
                size_factors=size_factor,
            )
        elif category == "igv":
            obj = ReadSegment.create(
                path=path,
                label=label,
                library=library,
                features=features,
                deletion_ignore=deletion_ignore,
                del_ratio_ignore=del_ratio_ignore,
                exon_focus=exon_focus,
            )
        elif category in ("hic", "h5"):
            obj = HiCTrack.create(
                path=path, label=label, log_trans=log_trans, depth=depth, tad=tad
            )
        elif category in ("bigwig", "bw"):
            category = "bw"
            obj = Bigwig.create(path, label=label, title=title)
        elif category in ("bedgraph", "bg"):
            category = "bg"
            obj = Bedgraph.create(path=path, label=label, title=title)
        elif category == "depth":
            obj = Depth.create(path, label=label, title=title)
        else:
            raise ValueError(
                f"the category should be one of [bam, bigwig, bw, depth, bedgraph, bg, h5], "
                f"instead of {category}"
            )
        obj.log_trans = log_trans
        return obj, category

    # ------------------------------------------------------------------
    # Public API: add different plot types
    # ------------------------------------------------------------------

    @add_object_error
    def add_customized_junctions(self, path: str):
        if path and os.path.exists(path) and os.path.isfile(path):
            self.junctions = load_custom_junction(path)

    @add_object_error
    def add_density(
        self,
        path,
        category="bam",
        label="",
        color="blue",
        font_size=8,
        show_junction_number=True,
        show_mean_jxn_number=False,
        junction_number_font_size=5,
        n_y_ticks=4,
        show_y_label=True,
        y_label="",
        theme="ticks_blank",
        show_site_plot=False,
        strand_choice=None,
        density_by_strand=False,
        only_customized_junction=False,
        **kwargs,
    ):
        obj, cat = self.__init_input_file__(
            path=path,
            category=category,
            label=label,
            density_by_strand=density_by_strand,
            **kwargs,
        )
        type_ = "site-plot" if show_site_plot and cat == "bam" else "density"
        if show_site_plot and cat != "bam":
            raise ValueError("show_site_plot only works with bam files")

        info = PlotInfo(obj=obj, type_=type_, category=cat)
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
                "show_mean_jxn_number": show_mean_jxn_number,
                "junction_number_font_size": junction_number_font_size,
                "color": color,
                "font_size": font_size,
                "n_y_ticks": n_y_ticks,
                "show_y_label": show_y_label,
                "y_label": y_label,
                "theme": theme,
                "strand_choice": strand_choice,
                "density_by_strand": density_by_strand,
                "only_customized_junction": only_customized_junction,
            }
        return self

    @add_object_error
    def add_heatmap(
        self,
        path,
        group="",
        category="bam",
        label="",
        color="viridis",
        font_size=8,
        show_y_label=True,
        theme="ticks_blank",
        do_scale=False,
        clustering=False,
        clustering_method="ward",
        distance_metric="euclidean",
        show_row_names=False,
        vmin=None,
        vmax=None,
        **kwargs,
    ):
        obj, cat = self.__init_input_file__(
            path=path, category=category, label=label, **kwargs
        )
        exists = False
        for p in self.plots:
            if p.group == group and p.type == "heatmap":
                p.add(obj=obj, category=cat, type_="heatmap")
                exists = True
                break
        if not exists:
            info = PlotInfo(obj=obj, category=cat, type_="heatmap", group=group)
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
                "vmin": vmin,
                "vmax": vmax,
            }
        return self

    @add_object_error
    def add_line(
        self,
        path,
        group="",
        category="bam",
        label="",
        color="blue",
        font_size=8,
        show_y_label=True,
        line_attrs=None,
        theme="ticks_blank",
        n_y_ticks=4,
        show_legend=False,
        legend_position="upper right",
        legend_ncol=0,
        **kwargs,
    ):
        obj, cat = self.__init_input_file__(
            path=path, category=category, label=label, **kwargs
        )
        if line_attrs is None:
            line_attrs = {}
        line_attrs["color"] = color

        exists = False
        for p in self.plots:
            if p.group == group and p.type == "line":
                params = self.params.pop(p)
                p.add(obj=obj, category=cat, type_="line")
                params["line_attrs"][obj.label] = line_attrs
                self.params[p] = params
                exists = True
                break
        if not exists:
            info = PlotInfo(obj=obj, category=cat, type_="line", group=group)
            self.plots.append(info)
            self.params[info] = {
                "line_attrs": {obj.label: line_attrs},
                "show_legend": show_legend,
                "font_size": font_size,
                "n_y_ticks": n_y_ticks,
                "show_y_label": show_y_label,
                "theme": theme,
                "legend_position": legend_position,
                "legend_ncol": legend_ncol,
            }
        return self

    @add_object_error
    def add_hic(
        self,
        path,
        category="hic",
        label="",
        color="RdYlBu_r",
        show_legend=True,
        font_size=8,
        n_y_ticks=4,
        show_y_label=True,
        theme="ticks",
        **kwargs,
    ):
        obj, cat = self.__init_input_file__(
            path=path, category=category, label=label, **kwargs
        )
        info = PlotInfo(obj=obj, category=cat, type_="hic")
        self.plots.append(info)
        self.params[info] = {
            "show_legend": show_legend,
            "color": color,
            "y_label": label,
            "font_size": font_size,
            "n_y_ticks": n_y_ticks,
            "show_y_label": show_y_label,
            "theme": theme,
        }
        return self

    @add_object_error
    def add_igv(
        self,
        path,
        category="igv",
        label="",
        library="fru",
        features=None,
        deletion_ignore=True,
        del_ratio_ignore=0.5,
        exon_focus=None,
        exon_color=None,
        intron_color=None,
        feature_color=None,
        exon_width=0.3,
        height_scale=None,
        font_size=8,
        n_y_ticks=1,
        show_y_label=True,
        theme="ticks_blank",
        **kwargs,
    ):
        obj, cat = self.__init_input_file__(
            path=path,
            category="igv",
            label=label,
            library=library,
            features=features,
            deletion_ignore=deletion_ignore,
            del_ratio_ignore=del_ratio_ignore,
            exon_focus=exon_focus,
        )
        info = PlotInfo(obj=obj, category=cat, type_="igv")
        info.height_scale = height_scale
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
            "theme": theme,
        }
        return self

    @add_object_error
    def add_motif(
        self,
        path,
        category="motif",
        motif_region=None,
        width=0.8,
        theme="blank",
        **kwargs,
    ):
        if motif_region.start < self.region.start:
            motif_region.start = self.region.start
        if motif_region.end > self.region.end:
            motif_region.end = self.region.end
        obj = Motif.create(path, motif_region)
        info = PlotInfo(obj=obj, category=category, type_="motif")
        self.plots.append(info)
        self.params[info] = {"width": width, "theme": theme}
        return self

    @add_object_error
    def add_manual(
        self,
        data,
        image_type="line",
        label="",
        group="",
        color="blue",
        font_size=8,
        n_y_ticks=1,
        show_y_label=True,
        theme="ticks_blank",
        **kwargs,
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
                "line_attrs": {obj.label: {"color": color}},
            }
            self.params[info].update(kwargs)
        return self

    def merge_by_cell(self):
        plots = {}
        for p in self.plots:
            assert isinstance(p, PlotInfo)
            if p.category in ["density", "line"]:
                label = p.obj[0].label
                if label in plots:
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
        return copy.deepcopy(self)

    # ------------------------------------------------------------------
    # Plot orchestration
    # ------------------------------------------------------------------

    def _load_plot_data(self, n_jobs, *args, **kwargs):
        """Phase 1: Load data for all plots (with optional multiprocessing)."""
        logger.info(f"load data of {len(self.plots)} plots")

        cmds = []
        for p in self.plots:
            assert isinstance(p, PlotInfo), f"unrecognized data type: {type(p)}"

            if self.__n_objs__ / len(self.plots) >= n_jobs > 1:
                if isinstance(p.obj[0].label, list):
                    juncs = {}
                    for i in p.obj[0].label:
                        juncs.update(self.junctions.get(i, {}))
                    p.load(self.region, junctions=juncs, *args, **kwargs)
                else:
                    p.load(
                        self.region,
                        n_jobs,
                        junctions=self.junctions.get(p.obj[0].label, {}),
                        *args,
                        **kwargs,
                    )
            elif n_jobs > 1:
                temp = copy.deepcopy(kwargs)
                temp["region"] = (
                    self.region if p.type != "motif" else self.params[p].get("region")
                )
                if isinstance(p.obj[0].label, list):
                    juncs = {}
                    for i in p.obj[0].label:
                        juncs.update(self.junctions.get(i, {}))
                    temp["junctions"] = juncs
                else:
                    temp["junctions"] = self.junctions.get(p.obj[0].label, {})
                cmds.append([p, args, temp])

        if cmds:
            with Pool(max(1, min(n_jobs, len(self.plots)))) as pool:
                self.plots = pool.map(__load__, cmds)

        for p in self.plots:
            if n_jobs <= 1:
                if isinstance(p.obj[0].label, list):
                    juncs = {}
                    for i in p.obj[0].label:
                        juncs.update(self.junctions.get(i, {}))
                    p.load(self.region, junctions=juncs, *args, **kwargs)
                else:
                    p.load(
                        self.region,
                        junctions=self.junctions.get(p.obj[0].label, {}),
                        *args,
                        **kwargs,
                    )

    def _build_layout(self, annotation_scale, sc_height_ratio, *args, **kwargs):
        """Phase 2: Compute figure layout dimensions."""
        igv_height_scale = kwargs.get("igv_height_scale")
        stroke_scale = kwargs.get("stroke_scale", 0.25)
        plots_n_rows, plots_n_cols = 1, 1
        height_ratio = []

        if self.annotation is not None:
            logger.info("load annotation")
            self.annotation.load(self.region, *args, **self.params["annotation"])
            plots_n_rows += self.annotation.len(scale=annotation_scale)
        if self.stroke:
            plots_n_rows += int(max(len(self.stroke) * stroke_scale, 1))
        if self.link:
            plots_n_rows += len(self.link)
        if self.sequence is not None:
            logger.info("load sequence")
            self.sequence.load(self.region)

        for p in self.plots:
            n_rows, n_height = p.len(
                annotation_scale,
                sc_height_ratio=sc_height_ratio,
                igv_scale=igv_height_scale,
            )
            plots_n_rows += n_rows
            height_ratio += n_height
            if p.type in ("heatmap", "hic"):
                plots_n_cols = 2

        plots_n_rows = int(plots_n_rows)
        height_ratio += [1] * (plots_n_rows - len(height_ratio))
        logger.debug(f"plots n_rows={plots_n_rows}; n_cols = {plots_n_cols}")
        logger.info("init graph_coords")
        return plots_n_rows, plots_n_cols, height_ratio

    def _create_figure_and_grid(
        self, plots_n_rows, plots_n_cols, height_ratio, **kwargs
    ):
        """Create matplotlib figure and GridSpec."""
        exon_scale = kwargs.get("exon_scale", 1)
        intron_scale = kwargs.get("intron_scale", 0.5)
        dpi = kwargs.get("dpi", 300)
        width = kwargs.get("width", 0)
        height = kwargs.get("height", 0)

        if plots_n_cols > 1 and intron_scale != 1:
            logger.debug("heatmap require intron_scale = 1")
            intron_scale = 1

        self.graph_coords = init_graph_coords(
            self.region, self.exons, exon_scale=exon_scale, intron_scale=intron_scale
        )

        fig = (
            plt.figure(figsize=[width, height * sum(height_ratio)], dpi=dpi)
            if width and height
            else plt.figure(dpi=dpi)
        )

        if plots_n_cols > 1:
            gs = gridspec.GridSpec(
                plots_n_rows,
                plots_n_cols,
                height_ratios=height_ratio,
                width_ratios=(0.99, 0.01),
                wspace=0.01,
                hspace=0.15,
            )
        else:
            gs = gridspec.GridSpec(
                plots_n_rows,
                plots_n_cols,
                height_ratios=height_ratio,
                wspace=0.7,
                hspace=0.15,
            )
        return fig, gs

    @staticmethod
    def _load_default_y(y_limit_path):
        """Load per-label y-limit settings from file."""
        default_y = {}
        if not y_limit_path or not os.path.exists(y_limit_path):
            if y_limit_path:
                logger.warning(f"{y_limit_path} not exists")
            return default_y
        logger.info("load y-limit settings from: " + y_limit_path)
        with open(y_limit_path) as r:
            for line in r:
                if line.startswith("#"):
                    continue
                line = line.split()
                if len(line) > 1:
                    try:
                        default_y[line[0]] = [float(x) for x in line[1:]]
                    except Exception as err:
                        logger.warning(f"The y limit of {line[0]} is invalid: {err}")
        return default_y

    def _precompute_plot_y_limits(
        self, p, max_used_y_dict, min_used_y_dict, base_max_dict, same_y_by_groups, default_y, **kwargs
    ):
        """Compute and store y-limits for all objects in a plot.

        Also collects base_max (raw data maximum) per object path for use as
        a global arc-height reference in --same-y mode.
        """
        if p.type not in ("density", "site-plot", "line"):
            return

        for obj in p.obj:
            _max, _min, _base = precompute_y_limits(
                obj,
                data=obj.data,
                region=self.region,
                graph_coords=self.graph_coords,
                **kwargs,
            )
            base_max_dict[obj.path] = max(_base, base_max_dict.get(obj.path, 0))

            if obj.label in same_y_by_groups:
                key = same_y_by_groups[obj.label]
                max_used_y_dict[key] = max(_max, max_used_y_dict.get(key, 0))
                min_used_y_dict[key] = min(
                    _min if obj.data.minus is not None else _min,
                    min_used_y_dict.get(key, 0)
                    if obj.data.minus is None
                    else min_used_y_dict.get(obj.path, 0),
                )
                base_max_dict[key] = max(_base, base_max_dict.get(key, 0))

            if isinstance(obj.data, dict):
                for key, readDepth in obj.data.items():
                    _max, _min, _base = precompute_y_limits(
                        obj,
                        data=readDepth,
                        graph_coords=self.graph_coords,
                        region=self.region,
                        **kwargs,
                    )
                    max_used_y_dict[key] = max(_max, max_used_y_dict.get(key, 0))
                    min_used_y_dict[key] = min(
                        _min if readDepth.minus is not None else _min,
                        min_used_y_dict.get(key, 0)
                        if readDepth.minus is None
                        else min_used_y_dict.get(obj.path, 0),
                    )
                    base_max_dict[key] = max(_base, base_max_dict.get(key, 0))

                continue

            max_used_y_dict[obj.path] = max(_max, max_used_y_dict.get(obj.path, 0))
            min_used_y_dict[obj.path] = min(
                _min if obj.data.minus is not None else _min,
                min_used_y_dict.get(obj.path, 0)
                if obj.data.minus is None
                else min_used_y_dict.get(obj.path, 0),
            )

            if obj.label in default_y:
                max_used_y_dict[obj.path] = default_y[obj.label][0]
                if len(default_y[obj.label]) > 1:
                    min_used_y_dict[obj.path] = default_y[obj.label][1]
                if same_y_by_groups and obj.label in same_y_by_groups:
                    key = same_y_by_groups[obj.label]
                    max_used_y_dict[key] = max_used_y_dict[obj.path]
                    min_used_y_dict[key] = min_used_y_dict[obj.path]

    def _resolve_plot_y_limits(
        self, p, max_used_y_val, min_used_y_val, base_max_val, same_y_by_groups, default_y, **kwargs
    ):
        """Determine y-limits for a single plot based on same-y / default-y settings.

        Returns (max_y, min_y, global_base_max).
        global_base_max is the max of all base_max values (for --same-y mode),
        used as the arc height reference so all panels share consistent arc heights.
        """
        same_y_sc = kwargs.get("same_y_sc")
        same_y = kwargs.get("same_y")

        if same_y_sc and p.obj[0].is_single_cell:
            _bm = base_max_val.get(p.obj[0].path) if base_max_val else None
            return max_used_y_val.get(p.obj[0].path), min_used_y_val.get(p.obj[0].path), _bm

        if same_y_by_groups and p.obj[0].label in same_y_by_groups:
            key = same_y_by_groups[p.obj[0].label]
            _bm = base_max_val.get(key) if base_max_val else None
            return max_used_y_val.get(key), min_used_y_val.get(key), _bm

        if same_y and max_used_y_val:
            _bm = max(base_max_val.values()) if base_max_val else None
            return max(max_used_y_val.values()), min(min_used_y_val.values()), _bm

        if default_y and p.type in ("density", "site-plot", "line"):
            for obj in p.obj:
                if obj.label in default_y:
                    max_y = default_y[obj.label][0]
                    min_y = (
                        default_y[obj.label][1]
                        if len(default_y[obj.label]) > 1
                        else None
                    )
                    if same_y_by_groups and obj.label in same_y_by_groups:
                        max_y = max_used_y_val.get(obj.path, max_y)
                        min_y = min_used_y_val.get(obj.path, min_y)
                    return max_y, min_y, None
        return None, None, None

    # ------------------------------------------------------------------
    # Main plot method
    # ------------------------------------------------------------------

    def plot(
        self,
        output=None,
        annotation_scale=0.25,
        stroke_scale=0.25,
        dpi=300,
        width=0,
        height=0,
        raster=False,
        return_image=None,
        sc_height_ratio=None,
        distance_between_label_axis=0.3,
        n_jobs=1,
        fill_step="post",
        *args,
        **kwargs,
    ):
        if sc_height_ratio is None:
            sc_height_ratio = {"density": 0.2, "heatmap": 0.2}
        assert self.region is not None, "please set the plotting region first."

        igv_height_scale = kwargs.get("igv_height_scale")

        # Extract / pop custom params before they flow into downstream **kwargs
        title = kwargs.pop("title", None)
        no_title = kwargs.pop("no_title", False)
        junctions_on_top = kwargs.pop("junctions_on_top", False)

        # ====== Phase 1: Load data ======
        self._load_plot_data(n_jobs, *args, **kwargs)

        # ====== Phase 2: Build layout ======
        plots_n_rows, plots_n_cols, height_ratio = self._build_layout(
            annotation_scale,
            sc_height_ratio,
            *args,
            stroke_scale=stroke_scale,
            **kwargs,
        )

        # ====== Phase 3: Create figure and grid ======
        fig, gs = self._create_figure_and_grid(
            plots_n_rows,
            plots_n_cols,
            height_ratio,
            dpi=dpi,
            width=width,
            height=height,
            **kwargs,
        )

        # ====== Phase 4: Precompute global y-limits ======
        max_used_y_val, min_used_y_val, base_max_val = {}, {}, {}
        same_y_by_groups = {}
        default_y = self._load_default_y(kwargs.get("y_limit"))

        if kwargs.get("same_y") or kwargs.get("same_y_sc"):
            logger.info("--same-y is enabled")
            if kwargs.get("same_y_groups") and os.path.exists(kwargs["same_y_groups"]):
                logger.info(
                    "Reading the same y by groups from: %s" % kwargs["same_y_groups"]
                )
                with open(kwargs["same_y_groups"]) as r:
                    for line in r:
                        if line.startswith("#"):
                            continue
                        line = line.split()
                        same_y_by_groups[line[0]] = (
                            f"samy_y_groups_of_{line[1]} {datetime.now()}"
                        )

            for p in self.plots:
                self._precompute_plot_y_limits(
                    p,
                    max_used_y_val,
                    min_used_y_val,
                    base_max_val,
                    same_y_by_groups,
                    default_y,
                    **kwargs,
                )

        # ====== Phase 5: Draw all plots ======
        curr_idx = 0
        for p in self.plots:
            if p.type == "igv":
                n_rows = p.len(annotation_scale, igv_scale=igv_height_scale)[0]
                ax_var = plt.subplot(gs[curr_idx : curr_idx + n_rows, 0])
            else:
                ax_var = plt.subplot(gs[curr_idx, 0])

            if curr_idx == 0 and not no_title:
                ax_var.set_title(title or str(self.region), loc="left")

            max_y_val_, min_y_val_, base_max_ = self._resolve_plot_y_limits(
                p,
                max_used_y_val,
                min_used_y_val,
                base_max_val,
                same_y_by_groups,
                default_y,
                same_y_sc=kwargs.get("same_y_sc"),
                same_y=kwargs.get("same_y"),
            )

            if max_y_val_ is not None:
                max_y_val_ *= 1.1

            if min_y_val_ is not None:
                min_y_val_ *= 1.1

            # For --same-y mode: use global base_max as arc height reference
            # so all panels draw junctions at a consistent scale.
            global_arc_ref = base_max_

            logger.info(
                f"plotting {p.type} at idx: {curr_idx} with height_ratio: {height_ratio[curr_idx]}"
            )

            # Dispatch
            if p.type == "density":
                if isinstance(p.obj[0], Depth):
                    for key, readDepth in p.obj[0].data.items():
                        temp_params = self.params.get(p, {}).copy()
                        if "y_label" not in temp_params:
                            temp_params["y_label"] = key
                        plot_density(
                            ax=ax_var,
                            data=readDepth,
                            region=self.region,
                            graph_coords=self.graph_coords,
                            max_used_y_val=max_y_val_,
                            min_used_y_val=min_y_val_,
                            distance_between_label_axis=distance_between_label_axis,
                            raster=raster,
                            fill_step=fill_step,
                            junctions_on_top=junctions_on_top,
                            global_arc_ref=global_arc_ref,
                            **temp_params,
                        )
                        curr_idx += 1
                        ax_var = plt.subplot(gs[curr_idx, 0])
                    curr_idx -= 1
                else:
                    plot_density(
                        ax=ax_var,
                        obj=p.obj[0],
                        graph_coords=self.graph_coords,
                        max_used_y_val=max_y_val_,
                        min_used_y_val=min_y_val_,
                        distance_between_label_axis=distance_between_label_axis,
                        raster=raster,
                        fill_step=fill_step,
                        junctions_on_top=junctions_on_top,
                        global_arc_ref=global_arc_ref,
                        **self.params.get(p, {}),
                    )
            elif p.type == "hic":
                plot_hic(
                    ax=ax_var,
                    cbar_ax=plt.subplot(gs[curr_idx, 1]),
                    obj=p.obj,
                    distance_between_label_axis=distance_between_label_axis,
                    raster=raster,
                    **self.params.get(p, {}),
                )
            elif p.type == "site-plot":
                plot_density(
                    ax=ax_var,
                    obj=p.obj[0],
                    graph_coords=self.graph_coords,
                    max_used_y_val=max_y_val_,
                    min_used_y_val=min_y_val_,
                    distance_between_label_axis=distance_between_label_axis,
                    raster=raster,
                    junctions_on_top=junctions_on_top,
                    global_arc_ref=global_arc_ref,
                    **self.params.get(p, {}),
                )
                curr_idx += 1
                plot_site_plot(
                    plt.subplot(gs[curr_idx, 0]),
                    p.obj[0],
                    graph_coords=self.graph_coords,
                    raster=raster,
                    distance_between_label_axis=distance_between_label_axis,
                    **self.params.get(p, {}),
                )
            elif p.type == "heatmap":
                plot_heatmap(
                    ax=ax_var,
                    cbar_ax=plt.subplot(gs[curr_idx, 1]),
                    data=p.data,
                    y_label=p.group,
                    graph_coords=self.graph_coords,
                    raster=raster,
                    distance_between_label_axis=distance_between_label_axis,
                    **self.params.get(p, {}),
                )
            elif p.type == "line":
                plot_line(
                    ax=ax_var,
                    data=p.data,
                    y_label=p.group,
                    max_used_y_val=max_y_val_,
                    min_used_y_val=min_y_val_,
                    graph_coords=self.graph_coords,
                    distance_between_label_axis=distance_between_label_axis,
                    **self.params.get(p, {}),
                )
            elif p.type == "igv":
                plot_igv_like(
                    ax=ax_var,
                    obj=p.data,
                    graph_coords=self.graph_coords,
                    raster=raster,
                    distance_between_label_axis=distance_between_label_axis,
                    **self.params.get(p, {}),
                )
            elif p.type == "motif":
                plot_motif(
                    ax=ax_var,
                    obj=p.obj[0],
                    graph_coords=self.graph_coords,
                    **self.params[p],
                )
            else:
                raise ValueError(f"unknown plot type {p.type}")

            set_indicator_lines(
                ax=ax_var, sites=self.sites, graph_coords=self.graph_coords
            )
            set_focus(ax=ax_var, focus=self.focus, graph_coords=self.graph_coords)

            curr_idx += (
                1
                if p.type != "igv"
                else p.len(annotation_scale, igv_scale=igv_height_scale)[0]
            )

        # ====== Links ======
        if self.link:
            for link in self.link:
                logger.info(f"plotting links at idx: {curr_idx}")
                plot_links(
                    ax=plt.subplot(gs[curr_idx : curr_idx + 1, 0]),
                    data=link,
                    graph_coords=self.graph_coords,
                )
                curr_idx += 1

        # ====== X-axis ticks ======
        logger.info(f"plotting x-axis ticks at idx: {curr_idx}")
        set_x_ticks(
            ax=plt.subplot(gs[curr_idx, 0]),
            region=self.region,
            graph_coords=self.graph_coords,
            sequence=self.sequence.data if self.sequence else None,
            **kwargs,
        )
        curr_idx += 1

        # ====== Annotation ======
        if self.annotation is not None:
            n_anno_rows = self.annotation.len(scale=annotation_scale)
            logger.info(f"plotting annotation at idx: {curr_idx}")
            ax_var = plt.subplot(gs[curr_idx : curr_idx + n_anno_rows, 0])
            plot_annotation(
                ax=ax_var,
                obj=self.annotation,
                graph_coords=self.graph_coords,
                plot_domain=self.annotation.add_domain,
                distance_between_label_axis=distance_between_label_axis,
                **self.params["annotation"],
            )
            set_indicator_lines(
                ax=ax_var, sites=self.sites, graph_coords=self.graph_coords
            )
            set_focus(ax=ax_var, focus=self.focus, graph_coords=self.graph_coords)
            curr_idx += n_anno_rows

        # ====== Stroke ======
        if self.stroke:
            logger.info(f"plotting stroke at idx: {curr_idx}")
            ax_var = plt.subplot(gs[curr_idx:plots_n_rows, 0])
            plot_stroke(
                ax=ax_var, data=self.stroke, graph_coords=self.graph_coords, **kwargs
            )

        # ====== Output ======
        if output:
            logger.info(f"saving fig into {output}")
            fig.savefig(output, transparent=False, bbox_inches="tight")
        elif return_image:
            buf = io.BytesIO()
            if return_image == "png":
                FigureCanvasAgg(fig).print_png(buf)
            elif return_image == "pdf":
                FigureCanvasPdf(fig).print_pdf(buf)
            return buf
        else:
            plt.show()

        plt.close()
        logger.info("Plot done")
        return self


if __name__ == "__main__":
    pass

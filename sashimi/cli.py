#!/usr/bin/env python3
# -*- coding: utf-8 -*-
u"""
Created at 2022.05.31

This script contains all the command line parameters
"""
import gzip
import os
import sys
from multiprocessing import cpu_count
from typing import Optional, Dict, Set, Tuple

import click
import matplotlib as mpl
from click_option_group import optgroup
from loguru import logger

from sashimi.base.GenomicLoci import GenomicLoci
from sashimi.conf.config import CLUSTERING_METHOD, COLORS, COLORMAP, DISTANCE_METRIC, IMAGE_TYPE
from sashimi.file.ATAC import ATAC
from sashimi.plot import Plot

__version__ = "0.1.2"
__author__ = "ygidtu & Ran Zhou"
__email__ = "ygidtu@gmail.com"


def decode_region(region: str):
    regions = region.split(":")

    if len(regions) < 3:
        strand = "+"
    else:
        strand = regions[-1]

    sites = [int(x) for x in regions[1].split("-")]
    return GenomicLoci(regions[0], sites[0], sites[1], strand)


class FileList(object):
    def __init__(self,
                 path: str,
                 category: str,
                 color: str = "black",
                 label: Optional[str] = None,
                 group: Optional[str] = None,
                 exon_focus: Optional[str] = None,
                 library: str = "fru",
                 trans: Optional[str] = None,
                 depth: int = 30000):

        self.path = os.path.abspath(os.path.expanduser(path))

        if not os.path.exists(self.path):
            raise FileNotFoundError(f"{self.path} not found.")

        self.label = label if label else os.path.basename(path)
        self.group = group
        self.color = color
        self.category = category
        self.exon_focus = exon_focus
        self.library = library
        self.trans = trans
        self.depth = depth

    @property
    def name(self) -> str:
        return os.path.basename(self.path)

    def __str__(self):
        return f"path: {self.path} \nlabel: {self.label} \ngroup: {self.group} \n" \
               f"color: {self.color} \ncategory: {self.category} \nlibrary: {self.library}"


def load_barcodes(barcode: str) -> Tuple[Dict[str, Dict[str, Set[str]]], Dict[str, str]]:
    u"""
    as name says
    :param barcode: the path to barcode file
    """
    r = gzip.open(barcode, "rt") if barcode.endswith(".gz") else open(barcode, "r")
    res = {}
    colors = {}
    for line in r:
        line = line.strip().split()
        if len(line) >= 4:
            key, bc, group, color = line[:4]
            colors[group] = color
        elif len(line) == 3:
            key, bc, group = line[:3]
        elif len(line) == 2:
            key, bc = line
            group = ""
        else:
            continue

        if key not in res.keys():
            res[key] = {}

        if group not in res[key]:
            res[key][group] = set()

        res[key][group].add(bc)

    r.close()
    return res, colors


def __read_iter__(path):
    with open(path) as r:
        for idx, line in enumerate(r):
            if line.startswith("#"):
                continue
            line = line.split()
            if not line:
                continue
            yield idx, line


def process_file_list(infile: str, category: str = "density"):
    u"""
    Process and check the file list format_
    :param infile: path to input file list
    :param category: the image type of file list used for
    """

    try:
        if category in ["density"]:
            for idx, line in __read_iter__(infile):
                path, category = line[0], line[1]

                if category not in ["bam", "bigwig", "bw", "depth", "igv", "atac"]:
                    raise ValueError(f"{category} is not supported in density plot.")

                if len(line) < 3:
                    yield FileList(path=path, category=category, color=COLORS[idx % len(COLORS)])
                elif len(line) < 4:
                    yield FileList(path=path, category=category, color=COLORS[idx % len(COLORS)], label=line[2])
                elif len(line) < 5:
                    yield FileList(path=path, category=category, color=line[3], label=line[2])
                elif len(line) < 6:
                    yield FileList(path=path, category=category, color=line[3], label=line[2], library=line[4])
                else:
                    yield FileList(path=path, category=category, color=line[3], label=line[2], library=line[4],
                                   depth=line[5])
        elif category in ["heatmap"]:
            groups = {}
            for idx, line in __read_iter__(infile):
                path, category = line[0], line[1]

                if category not in ["bam", "bigwig", "bw", "depth", "atac"]:
                    raise ValueError(f"{category} is not supported in heatmap plot.")

                if len(line) < 3:
                    yield FileList(path=path, category=category, color=COLORMAP[0])
                elif len(line) < 4:
                    groups[line[2]] = 0
                    yield FileList(path=path, category=category,
                                   color=COLORMAP[len(groups) % len(COLORMAP)], group=line[2])
                else:
                    groups[line[2]] = 0
                    yield FileList(path=path, category=category,
                                   color=line[3], group=line[2], library=line[4] if len(line) > 4 else "fru")
        elif category in ["line"]:
            groups = {}
            for idx, line in __read_iter__(infile):
                path, category = line[0], line[1]

                if category not in ["bam", "bigwig", "bw", "depth"]:
                    raise ValueError(f"{category} is not supported in density plot.")

                if len(line) < 3:
                    yield FileList(path=path, category=category, color=COLORS[idx % len(COLORS)])
                elif len(line) < 4:
                    if line[2] not in groups:
                        groups[line[2]] = 0
                    groups[line[2]] += 1
                    yield FileList(path=path, category=category,
                                   color=COLORS[groups[line[2]] % len(COLORS)], group=line[2])
                elif len(line) < 5:
                    if line[2] not in groups:
                        groups[line[2]] = 0
                    groups[line[2]] += 1
                    yield FileList(path=path, category=category, label=line[3],
                                   color=COLORS[groups[line[2]] % len(COLORS)], group=line[2])
                elif len(line) < 6:
                    if line[2] not in groups:
                        groups[line[2]] = 0
                    groups[line[2]] += 1
                    yield FileList(path=path, category=category, label=line[3],
                                   color=line[4], group=line[2])
        elif category in ["interval"]:
            for idx, line in __read_iter__(infile):
                if len(line) < 2:
                    yield FileList(path=line[0], category="interval")
                else:
                    yield FileList(path=line[0], category="interval", label=line[1])
        elif category in ["igv"]:
            for idx, line in __read_iter__(infile):
                path, category = line[0], line[1]

                if category not in ["bam", "bigwig", "bw", "depth", "igv"]:
                    raise ValueError(f"{category} is not supported in density plot.")

                if len(line) < 3:
                    yield FileList(path=path, category=category, color=COLORS[idx % len(COLORS)])
                elif len(line) < 4:
                    yield FileList(path=path, category=category, color=COLORS[idx % len(COLORS)], label=line[2])
                elif len(line) < 5:
                    yield FileList(path=path, category=category, color=line[3], label=line[2])
                else:
                    yield FileList(path=path, category=category, color=line[3], label=line[2], exon_focus=line[4])
        elif category in ["hic"]:
            for idx, line in __read_iter__(infile):
                path, category = line[0], line[1]
                if len(line) < 3:
                    yield FileList(path=path, category=category)
                elif len(line) < 4:
                    yield FileList(path=path, category=category, label=line[2])
                elif len(line) < 5:
                    yield FileList(path=path, category=category, )
                else:
                    yield FileList(path=path, category=category, label=line[2], color=line[3], trans=line[4])
    except FileNotFoundError as err:
        logger.error(f"{infile} -> {err}")
        exit(1)

    return None


@click.command(context_settings=dict(help_option_names=['-h', '--help']), )
@click.version_option(__version__, message="Current version %(version)s")
@click.option("--debug", is_flag=True, help="enable debug level log")
@click.option("-e", "--event", type=click.STRING, required=True,
              help="Event range eg: chr1:100-200:+")
@optgroup.group("Common input files configuration")
@optgroup.option("--color-factor", default=1, type=click.IntRange(min=1),
                 help="Index of column with color levels (1-based); "
                      "NOTE: LUAD|red -> LUAD while be labeled in plots and red while be the fill color",
                 show_default=True)
@optgroup.option("--barcode", type=click.Path(exists=True), show_default=True,
                 help="Path to barcode list file, At list two columns were required, "
                      "- 1st The name of bam file, not the alias of bam; \b"
                      "- 2nd the barcode; \b"
                      "- 3rd The group label, optional; \b"
                      "- 4th The color of each cell type, default using the color of corresponding bam file.\n")
@optgroup.option("--barcode-tag", type=click.STRING, default="CB", show_default=True,
                 help="The default cell barcode tag label")
@optgroup.option("--umi-tag", type=click.STRING, default="UB", show_default=True,
                 help="The default UMI barcode tag label")
@optgroup.option("-p", "--process", type=click.IntRange(min=1, max=cpu_count()), default=1,
                 help="How many cpu to use")
@optgroup.option("--group-by-cell", type=click.BOOL, is_flag=True, help="Group by cell types in density/line plot")
@optgroup.option("--remove-duplicate-umi", type=click.BOOL, is_flag=True, help="Drop duplicated UMIs by barcode")
@optgroup.group("Output settings")
@optgroup.option("-o", "--output", type=click.Path(),
                 help="Path to output graph file", show_default=True)
@optgroup.option("-d", "--dpi", default=300, type=click.IntRange(min=1, clamp=True),
                 help="The resolution of output file", show_default=True)
@optgroup.option("--raster", is_flag=True, show_default=True,
                 help="The would convert heatmap and site plot to raster image "
                      "(speed up rendering and produce smaller files), only affects pdf, svg and PS")
@optgroup.option("--height", default=1, type=float,
                 help="The height of output file, default adjust image height by content", show_default=True)
@optgroup.option("--width", default=10, type=click.IntRange(min=0, clamp=True),
                 help="The width of output file, default adjust image width by content", show_default=True)
@optgroup.option("--backend", type=click.STRING, default="Agg", help="Recommended backend", show_default=True)
@optgroup.group("Reference settings")
@optgroup.option("-r", "--reference", type=click.Path(exists=True),
                 help="Path to gtf file, both transcript and exon tags are necessary")
@optgroup.option("--interval", type=click.Path(exists=True),
                 help="Path to list of interval files in bed format, "
                      "1st column is path to file, 2nd column is the label [optional]")
@optgroup.option("--show-id", is_flag=True, show_default=True, help="Whether show gene id or gene name")
@optgroup.option("--show-exon-id", is_flag=True, show_default=True, help="Whether show gene id or gene name")
@optgroup.option("--no-gene", is_flag=True, type=click.BOOL, show_default=True,
                 help="Do not show gene id next to transcript id")
@optgroup.option("--domain", default=False, is_flag=True, type=click.BOOL, show_default=True,
                 help="Add domain information into reference track")
@optgroup.option("--proxy", default=None, type=click.STRING,
                 help="The http or https proxy for EBI/Uniprot requests,"
                      "if `--domain` is True, eg: http://127.0.0.1:1080")
@optgroup.option("--timeout", default=10, type=click.IntRange(min=1, clamp=True),
                 show_default=True,
                 help="The requests timeout when `--domain` is True.")
@optgroup.option("--local-domain", default="", is_flag=False, type=click.STRING, show_default=True,
                 help="Load local domain folder and load into reference track, download from "
                      "https://hgdownload.soe.ucsc.edu/gbdb/hg38/uniprot/")
@optgroup.option("--remove-empty", is_flag=True, type=click.BOOL, show_default=True,
                 help="Whether to plot empty transcript")
@optgroup.option("--transcripts-to-show", default="", show_default=True,
                 help="Which transcript to show, transcript name or id in gtf file, eg: transcript1,transcript2")
@optgroup.option("--choose-primary", is_flag=True, type=click.BOOL, show_default=True,
                 help="Whether choose primary transcript to plot.")
@optgroup.option("--ref-color", default="black", type=click.STRING,
                 show_default=True, help="The color of exons")
@optgroup.option("--intron-scale", type=click.FLOAT, default=0.5, help="The scale of intron", show_default=True)
@optgroup.option("--exon-scale", type=click.FLOAT, default=1, help="The scale of exon", show_default=True)
@optgroup.group("Density plot settings")
@optgroup.option("--density", type=click.Path(exists=True),
                 help="""
                 The path to list of input files, a tab separated text file, \b
                 - 1st column is path to input file, \b
                 - 2nd column is the file category, \b
                 - 3rd column is input file alias (optional), \b
                 - 4th column is color of input files (optional),
                 - 5th column is the library of input file (optional, only required by bam file). \n
                 """)
@optgroup.option("--customized-junction", type=click.STRING, default=None, show_default=True,
                 help="Path to junction table column name needs to be bam name or bam alias.")
@optgroup.option("--only-customized-junction", is_flag=True, show_default=True, help="Only used customized junctions.")
@optgroup.option("-t", "--threshold", default=0, type=click.IntRange(min=0, clamp=True),
                 show_default=True, help="Threshold to filter low abundance junctions")
@optgroup.option("--density-by-strand", is_flag=True, type=click.BOOL,
                 show_default=True, help="Whether to draw density plot by strand")
@optgroup.option("--show-site", is_flag=True, type=click.BOOL,
                 show_default=True, help="Whether to draw additional site plot")
@optgroup.option("--site-strand", type=click.Choice(["all", "+", "-"]), default="all", show_default=True,
                 help="Which strand kept for site plot, default use all")
@optgroup.option("--included-junctions", type=click.STRING, default=None,
                 help="The junction id for including, chr1:1-100", show_default=True)
@optgroup.option("--show-junction-num", type=click.BOOL, is_flag=True, show_default=True,
                 help="Whether to show the number of junctions")
@optgroup.option("--sc-density-height-ratio", type=float, default=1, show_default=True,
                 help="The relative height of single cell density plots")
@optgroup.group("Line plot settings")
@optgroup.option("--line", type=click.Path(exists=True),
                 help="""
                 The path to list of input files, a tab separated text file, \b
                 - 1st column is path to input file, \b
                 - 2nd column is the file category, \b
                 - 3rd column is input file group (optional), \b
                 - 4th column is input file alias (optional),
                 - 5th column is color platte of corresponding group (optional).
                 """)
@optgroup.option("--hide-legend", default=False, is_flag=True, type=click.BOOL, help="Whether to hide legend")
@optgroup.option("--legend-position", default="upper right", type=click.STRING, help="The legend position")
@optgroup.option("--legend-ncol", default=0, type=click.IntRange(min=0, clamp=True),
                 help="The number of columns of legend")
@optgroup.group("Heatmap plot settings")
@optgroup.option("--heatmap", type=click.Path(exists=True),
                 help="""
                 The path to list of input files, a tab separated text file, \b
                 - 1st column is path to input file, \b
                 - 2nd column is the file category, \b
                 - 3rd column is input file group (optional), \b
                 - 4th column is color platte of corresponding group.
                 """)
@optgroup.option("--clustering", is_flag=True, show_default=True, help="Enable clustering of the heatmap")
@optgroup.option("--clustering-method", type=click.Choice(CLUSTERING_METHOD), default="ward",
                 show_default=True, help="The clustering method for heatmap")
@optgroup.option("--distance-metric", type=click.Choice(DISTANCE_METRIC), default="euclidean",
                 show_default=True, help="The distance metric for heatmap")
@optgroup.option("--heatmap-scale", is_flag=True, show_default=True, help="Do scale on heatmap matrix.")
@optgroup.option("--heatmap-vmin", type=click.INT, show_default=True,
                 help="Minimum value to anchor the colormap, otherwise they are inferred from the data.")
@optgroup.option("--heatmap-vmax", type=click.INT, show_default=True,
                 help="Maximum value to anchor the colormap, otherwise they are inferred from the data.")
@optgroup.option("--show-row-names", is_flag=True, show_default=True, help="Show row names of heatmap")
@optgroup.option("--sc-heatmap-height-ratio", type=float, default=.2, show_default=True,
                 help="The relative height of single cell heatmap plots")
@optgroup.group("IGV settings")
@optgroup.option("--igv", type=click.Path(exists=True),
                 help="""
                 The path to list of input files, a tab separated text file, \b
                 - 1st column is path to input file, \b
                 - 2nd column is the file category, \b
                 - 3rd column is input file alias (optional), \b
                 - 4th column is color of input files (optional),\b
                 - 5th column is exon_id for sorting the reads (optional).
                 """)
@optgroup.option("--m6a", default=None, type=click.STRING,
                 help="""
                 Sashimi.py will load location information from the given tags and \b
                  then highlight the RNA m6a modification cite at individual reads. \b
                  If there are multiple m6a modification site, please add tag as follow, \b
                  234423,234450
                 """)
@optgroup.option("--polya", default=None, type=click.STRING,
                 help="""
                 Sashimi.py will load length of poly(A) from the given tags and \b
                  then visualize the poly(A) part at end of each individual reads.
                 """)
@optgroup.option("--rs", default=None, type=click.STRING,
                 help="""
                 Sashimi.py will load real strand information of each reads from the given tags and \b
                  the strand information is necessary for visualizing poly(A) part.
                 """)
@optgroup.option("--del-ratio-ignore", default=1.0,
                 type=click.FloatRange(min=0.0, max=1.0, clamp=True),
                 help="""
                 Ignore the deletion gap in nanopore or pacbio reads. \b
                 if a deletion region was smaller than (alignment length) * (del_ratio_ignore), \b
                 then the deletion gap will be filled. \b
                 currently the del_ratio_ignore was 1.0.
                 """)
@optgroup.group("HiC settings")
@optgroup.option("--hic", type=click.Path(exists=True),
                 help="""
                 The path to list of input files, a tab separated text file, \b
                 - 1st column is path to input file, \b
                 - 2nd column is the file category, \b
                 - 3rd column is input file alias (optional), \b
                 - 4th column is color of input files (optional),\b
                 - 5th column is data transform for HiC matrix, eg log1p, log2, log10 (optional).
                 """)
@optgroup.group("Additional annotation")
@optgroup.option("-f", "--genome", type=click.Path(), default=None,
                 show_default=True, help="Path to genome fasta")
@optgroup.option("--sites", default=None, type=click.STRING,
                 help="Where to plot additional indicator lines, comma separated int")
@optgroup.option("--stroke", type=click.STRING, show_default=True,
                 help="The stroke regions: start1-end1:start2-end2@color-label, "
                      "draw a stroke line at bottom, default color is red")
@optgroup.option("--link", type=click.STRING, show_default=True,
                 help="The link: start1-end1:start2-end2@color, "
                      "draw a link between two site at bottom, default color is blue")
@optgroup.option("--focus", type=click.STRING, show_default=True, help="The highlight regions: 100-200:300-400")
@optgroup.group("Motif settings")
@optgroup.option("--motif", type=click.Path(exists=True),
                 help="The path to customized bedGraph file, first three columns is chrom, start and end site, "
                      "the following 4 columns is the weight of ATCG.")
@optgroup.option("--motif-region", type=click.STRING, default="",
                 help="The region of motif to plot in start-end format", show_default=True)
@optgroup.option("--motif-width", type=click.FLOAT, default=0.8,
                 help="The width of ATCG characters", show_default=True)
@optgroup.group("Layout settings")
@optgroup.option("--n-y-ticks", default=4, type=click.IntRange(min=0, clamp=True),
                 help="The number of ticks of y-axis")
@optgroup.option("--distance-ratio", type=click.FLOAT, default=0.1,
                 help="distance between transcript label and transcript line", show_default=True)
@optgroup.option("--reference-scale", type=click.FLOAT, default=.25,
                 help="The size of reference plot in final plot", show_default=True)
@optgroup.option("--stroke-scale", type=click.FLOAT, default=.25,
                 help="The size of stroke plot in final image", show_default=True)
@optgroup.group("Overall settings")
@optgroup.option("--font-size", default=8, type=click.IntRange(min=1, clamp=True),
                 help="The font size of x, y-axis and so on")
@optgroup.option("--reverse-minus", default=False, is_flag=True, type=click.BOOL,
                 help="Whether to reverse strand of bam/reference file")
@optgroup.option("--hide-y-label", default=False, is_flag=True, type=click.BOOL,
                 help="Whether hide y-axis label")
@optgroup.option("--same-y", default=False, is_flag=True, type=click.BOOL,
                 help="Whether different sashimi/line plots shared same y-axis boundaries")
@optgroup.option('--log', type=click.Choice(["0", "2", "10", "zscore"]), default="0",
                 help="y axis log transformed, 0 -> not log transform; 2 -> log2; 10 -> log10")
@optgroup.option("--title", type=click.STRING, default=None, help="Title", show_default=True)
@optgroup.option("--font", type=click.STRING, default=None, help="Fonts", show_default=True)
def main(**kwargs):
    u"""
    Welcome to use sashimi
    \f
    """
    if not kwargs["debug"]:
        logger.remove()
        logger.add(sys.stderr, level="INFO")
        logger.level("INFO")

    # print warning info about backend
    if (kwargs["domain"] or kwargs["local_domain"]) and kwargs["backend"].lower() != "cairo":
        logger.warning(f"{kwargs['backend']} backend may have problems with small domain, "
                       f"if there is any please try cairo backend instead.")

    if kwargs["raster"] and kwargs["heatmap"] and kwargs["backend"].lower() == "cairo":
        logger.warning(f"{kwargs['backend']} backend may have problems with rasterized heatmap, "
                       f"if there is any, please try another backend instead.")

    try:
        mpl.use(kwargs["backend"])
    except ImportError as err:
        if kwargs["backend"].lower() == "cairo":
            logger.warning("Cairo backend required cairocffi installed")
            logger.warning("Switch back to Agg backend")
        else:
            logger.warning(f"backend error, switch back to Agg: {err}")
        mpl.use("Agg")

    mpl.rcParams['pdf.fonttype'] = 42

    if kwargs["font"]:
        mpl.rcParams['font.family'] = kwargs["font"]

    for k, v in kwargs.items():
        logger.debug(f"{k} => {v}")

    if kwargs["included_junctions"] is not None:
        included_junctions = set([sj.strip() for sj in kwargs["included_junctions"].split(',')])
    else:
        included_junctions = {}

    p = Plot()

    region = decode_region(kwargs["event"])
    p.set_region(region=region)
    p.add_customized_junctions(kwargs["customized_junction"])

    barcodes, sc_colors = {}, {}
    if kwargs.get("barcode") and os.path.exists(kwargs.get("barcode")):
        barcodes, sc_colors = load_barcodes(kwargs.get("barcode"))

    size_factors = {}

    # add reference
    for key in kwargs.keys():
        # print(key)
        if key in IMAGE_TYPE and kwargs[key] and os.path.exists(kwargs[key]):
            logger.debug(f"add {key} {kwargs[key]}")
            if key == "reference":
                p.set_reference(kwargs["reference"],
                                show_gene=not kwargs["no_gene"],
                                color=kwargs["ref_color"],
                                remove_empty_transcripts=kwargs["remove_empty"],
                                choose_primary=kwargs["choose_primary"],
                                font_size=kwargs["font_size"],
                                show_id=kwargs["show_id"],
                                reverse_minus=kwargs["reverse_minus"],
                                show_exon_id=kwargs["show_exon_id"],
                                transcripts=kwargs["transcripts_to_show"],
                                add_domain=kwargs["domain"],
                                local_domain=kwargs["local_domain"]
                                )
            elif key == "interval":
                for f in process_file_list(kwargs[key], key):
                    p.add_interval(f.path, f.label)
            elif key == "density":
                for f in process_file_list(kwargs[key], key):
                    if barcodes and f.name in barcodes.keys() and f.category in ["bam", "atac"]:
                        for group in barcodes[f.name].keys():
                            if kwargs["group_by_cell"] and group:
                                label = group
                            elif group:
                                label = f"{f.label} - {group}"
                            else:
                                label = f.label

                            if f.label not in size_factors.keys() and f.category == "atac":
                                logger.info(f"Indexing {f.path}")
                                size_factors[f.label] = ATAC.index(f.path, barcodes[f.name])

                            p.add_density(f.path,
                                          category=f.category,
                                          label=label,
                                          barcode=group,
                                          barcode_groups=barcodes[f.name],
                                          barcode_tag=kwargs["barcode_tag"],
                                          umi_tag=kwargs["umi_tag"],
                                          library=f.library,
                                          size_factor=size_factors.get(f.label),
                                          color=sc_colors.get(group, f.color),
                                          font_size=kwargs["font_size"],
                                          show_junction_number=kwargs["show_junction_num"],
                                          n_y_ticks=kwargs["n_y_ticks"],
                                          show_y_label=not kwargs["hide_y_label"],
                                          show_site_plot=kwargs["show_site"],
                                          strand_choice=kwargs["site_strand"],
                                          density_by_strand=kwargs["density_by_strand"],
                                          only_customized_junction=kwargs["only_customized_junction"],
                                          log_trans=kwargs["log"])
                    elif f.category != "atac":
                        p.add_density(f.path,
                                      category=f.category,
                                      label=f.label,
                                      barcode_tag=kwargs["barcode_tag"],
                                      umi_tag=kwargs["umi_tag"],
                                      library=f.library,
                                      color=f.color,
                                      font_size=kwargs["font_size"],
                                      show_junction_number=kwargs["show_junction_num"],
                                      n_y_ticks=kwargs["n_y_ticks"],
                                      show_y_label=not kwargs["hide_y_label"],
                                      show_site_plot=kwargs["show_site"],
                                      strand_choice=kwargs["site_strand"],
                                      density_by_strand=kwargs["density_by_strand"],
                                      log_trans=kwargs["log"])
            elif key == "heatmap":
                for f in process_file_list(kwargs[key], key):
                    if barcodes and f.name in barcodes.keys() and f.category in ["bam", "atac"]:
                        if f.label not in size_factors.keys() and f.category == "atac":
                            logger.info(f"Indexing {f.path}")
                            size_factors[f.label] = ATAC.index(f.path, barcodes[f.name])

                        for group in barcodes[f.name].keys():
                            p.add_heatmap(f.path,
                                          category=f.category,
                                          label=f"{f.label} - {group}" if group else f.label,
                                          barcode=group,
                                          barcode_groups=barcodes[f.name],
                                          group=f"{f.group} - {group}" if f.group else f.group,
                                          barcode_tag=kwargs["barcode_tag"],
                                          size_factor=size_factors.get(f.label),
                                          umi_tag=kwargs["umi_tag"],
                                          library=f.library,
                                          color=f.color,
                                          show_y_label=not kwargs["hide_y_label"],
                                          clustering=kwargs["clustering"],
                                          clustering_method=kwargs["clustering_method"],
                                          distance_metric=kwargs["distance_metric"],
                                          font_size=kwargs["font_size"],
                                          do_scale=kwargs["heatmap_scale"],
                                          vmin=kwargs["heatmap_vmin"],
                                          vmax=kwargs["heatmap_vmax"],
                                          log_trans=kwargs["log"])
                    elif f.category != "atac":
                        p.add_heatmap(f.path,
                                      category=f.category,
                                      group=f.group,
                                      label=f.label,
                                      barcode_tag=kwargs["barcode_tag"],
                                      umi_tag=kwargs["umi_tag"],
                                      library=f.library,
                                      color=f.color,
                                      show_y_label=not kwargs["hide_y_label"],
                                      clustering=kwargs["clustering"],
                                      clustering_method=kwargs["clustering_method"],
                                      distance_metric=kwargs["distance_metric"],
                                      font_size=kwargs["font_size"],
                                      show_row_names=kwargs["show_row_names"],
                                      do_scale=kwargs["heatmap_scale"],
                                      vmin=kwargs["heatmap_vmin"],
                                      vmax=kwargs["heatmap_vmax"],
                                      log_trans=kwargs["log"])
            elif key == "line":
                for f in process_file_list(kwargs[key], key):
                    if barcodes and f.name in barcodes.keys() and f.category == "bam":
                        for group in barcodes[f.name].keys():
                            if kwargs["group_by_cell"] and group:
                                label = group
                            elif group:
                                label = f"{f.label} - {group}"
                            else:
                                label = f.label
                            p.add_line(f.path,
                                       category=f.category,
                                       label=label,
                                       barcode=group,
                                       barcode_groups=barcodes,
                                       group=f.group,
                                       barcode_tag=kwargs["barcode_tag"],
                                       umi_tag=kwargs["umi_tag"],
                                       library=f.library,
                                       color=sc_colors.get(group, f.color),
                                       show_y_label=not kwargs["hide_y_label"],
                                       font_size=kwargs["font_size"],
                                       n_y_ticks=kwargs["n_y_ticks"],
                                       show_legend=not kwargs["hide_legend"],
                                       legend_position=kwargs["legend_position"],
                                       legend_ncol=kwargs["legend_ncol"],
                                       log_trans=kwargs["log"])
                    else:
                        p.add_line(f.path,
                                   category=f.category,
                                   group=f.group,
                                   label=f.label,
                                   barcode_tag=kwargs["barcode_tag"],
                                   umi_tag=kwargs["umi_tag"],
                                   library=f.library,
                                   color=f.color,
                                   show_y_label=not kwargs["hide_y_label"],
                                   font_size=kwargs["font_size"],
                                   n_y_ticks=kwargs["n_y_ticks"],
                                   show_legend=not kwargs["hide_legend"],
                                   legend_position=kwargs["legend_position"],
                                   legend_ncol=kwargs["legend_ncol"],
                                   log_trans=kwargs["log"])
            elif key == "igv":
                for f in process_file_list(kwargs[key], "igv"):
                    igv_features = {}
                    if kwargs["m6a"]:
                        igv_features.update({"m6a": kwargs["m6a"]})

                    if kwargs["polya"] and kwargs["rs"]:
                        igv_features.update({"real_strand": kwargs["rs"], "polya": kwargs["polya"]})

                    if len(igv_features) == 0:
                        igv_features = None

                    p.add_igv(f.path,
                              category=f.category,
                              label=f.label,
                              exon_color=f.color,
                              intron_color=f.color,
                              features=igv_features,
                              font_size=kwargs["font_size"],
                              n_y_ticks=kwargs["n_y_ticks"],
                              show_y_label=not kwargs["hide_y_label"],
                              deletion_ignore=True if kwargs["del_ratio_ignore"] == 1.0 else False,
                              del_ratio_ignore=kwargs["del_ratio_ignore"],
                              exon_focus=f.exon_focus
                              )
            elif key == "hic":
                for f in process_file_list(kwargs[key], "hic"):
                    p.add_hic(
                        f.path,
                        category=f.category,
                        label=f.label,
                        trans=f.trans,
                        depth=f.depth,
                        color=f.color,
                        show_legend=not kwargs["hide_legend"],
                        show_y_label=not kwargs["hide_y_label"],
                        font_size=kwargs["font_size"],
                        n_y_ticks=kwargs["n_y_ticks"]
                    )
            elif key == "motif":
                motif_region = None
                if kwargs["motif_region"]:
                    start, end = [int(x) for x in kwargs["motif_region"].split("-")]
                    motif_region = GenomicLoci(
                        region.chromosome,
                        max(start, region.start),
                        min(region.end, end),
                        region.strand)

                p.add_motif(kwargs[key], motif_region=motif_region, width=kwargs["motif_width"])
        elif key == "focus":
            p.add_focus(kwargs[key])
        elif key == "stroke":
            p.add_stroke(kwargs[key])
        elif key == "sites":
            p.add_sites(kwargs[key])
        elif key == "link":
            p.add_links(kwargs[key])

    if kwargs["group_by_cell"]:
        p.merge_by_cell()

    p.plot(
        kwargs["output"],
        width=kwargs["width"],
        height=kwargs["height"],
        dpi=kwargs["dpi"],
        raster=kwargs["raster"],
        intron_scale=kwargs["intron_scale"],
        exon_scale=kwargs["exon_scale"],
        reference_scale=kwargs["reference_scale"],
        strock_scale=kwargs["stroke_scale"],
        same_y=kwargs["same_y"],
        remove_duplicate_umi=kwargs["remove_duplicate_umi"],
        threshold=kwargs["threshold"],
        sc_height_ratio={
            "heatmap": kwargs["sc_heatmap_height_ratio"],
            "density": kwargs["sc_density_height_ratio"]
        },
        distance_between_label_axis=kwargs["distance_ratio"],
        included_junctions=included_junctions,
        n_jobs=kwargs.get("process", 1)
    )
    logger.info("DONE")


if __name__ == '__main__':
    main()

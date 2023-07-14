#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
Web UI of sashimi
"""
import pickle
import re
from glob import glob

import click

from flask import Flask, render_template, jsonify, send_file, request

from trackplot.cli import load_barcodes, __version__
from trackplot.plot import *

__DIR__ = os.path.abspath(os.path.dirname(__file__))
__UI__ = os.path.join(__DIR__, "ui")
__PLOT__ = os.path.join(os.path.dirname(__file__), "plots")

app = Flask( __name__, static_url_path="/static", static_folder=__UI__, template_folder=__UI__)


# supported trackplot settings
__COMMON_PARAMS__ = [
    {
        "key": "path",
        "annotation": "str",
        "default": "<class 'inspect._empty'>",
        "note": "Please input the file path"
    },
    {
        "key": "category",
        "annotation": "str",
        "default": "bam",
        "note": "The category of input file"
    },
    {
        "key": "size_factor",
        "annotation": "empty",
        "default": "<class 'inspect._empty'>",
        "note": "Total number of reads for bam file or total number of fragments required by scATAC"
    },
    {
        "key": "label",
        "annotation": "Union[str, List[str]]",
        "default": "<class 'inspect._empty'>",
        "note": "The alias of input file"
    },
    {
        "key": "title",
        "annotation": "str",
        "default": "<class 'inspect._empty'>",
        "note": "The title of this sub plot"
    },
    {
        "key": "barcode_groups",
        "annotation": "str",
        "default": "<class 'inspect._empty'>",
        "note": "The path to list of barcodes"
    },
    {
        "key": "library",
        "annotation": "choice['frr', 'frs', 'fru']",
        "default": "fru",
        "note": "The strand preference of input file, frf => fr-firststrand; frs => fr-secondstrand; fru => fr-unstrand"
    },
    {
        "key": "n_y_ticks",
        "annotation": "int",
        "default": "4",
        "note": "Number of y ticks"
    },
    {
        "key": "show_y_label",
        "annotation": "bool",
        "default": "true",
        "note": "Whether to show y label"
    },
    {
        "key": "y_label",
        "annotation": "str",
        "default": "<class 'inspect._empty'>",
        "note": "The title of y-axis"
    },
    {
        "key": "color",
        "annotation": "color",
        "default": "#409EFF",
        "note": "The fill color of density"
    },
    {
        "key": "font_size",
        "annotation": "int",
        "default": "8",
        "note": "The font size of this sub plot"
    },
    {
        "key": "log_trans",
        "annotation": "choice['0', '2', '10']",
        "default": "0",
        "note": "Whether to perform log transformation, 0 -> not log transform;2 -> log2;10 -> log10"
    },
]

__PARAMS__ = {
    "add_customized_junctions": [
        {
            "key": "path",
            "annotation": "str",
            "default": "<class 'inspect._empty'>",
            "note": "Path to junction table column name needs to be bam name or bam alias."
        }
    ],
    "add_density": [
        [
            "path", "category", "size_factor", "label", "title",
            "barcode_groups", "barcode_tag", "umi_tag", "library",
            "color", "font_size", "n_y_ticks", "show_y_label", "y_label",
            "log_trans"
        ],
        {
            "key": "density_by_strand",
            "annotation": "bool",
            "default": "false",
            "note": "Whether draw density plot by strand"
        },
        {
            "key": "show_junction_number",
            "annotation": "bool",
            "default": "true",
            "note": "Whether to show junction number"
        },
        {
            "key": "junction_number_font_size",
            "annotation": "int",
            "default": "5",
            "note": "The font size of junction number"
        },
        {
            "key": "show_site_plot",
            "annotation": "bool",
            "default": "false",
            "note": "Whether to show site plot"
        },
        {
            "key": "strand_choice",
            "annotation": "choice['all', '+', '-']",
            "default": "all",
            "note": "Which strand kept for site plot, default use all"
        },
        {
            "key": "only_customized_junction",
            "annotation": "bool",
            "default": "false",
            "note": "Only used customized junctions."
        }
    ],
    "add_focus": [
        {
            "key": "start",
            "annotation": "int",
            "default": "0",
            "note": "The start site of focus region"
        },
        {
            "key": "end",
            "annotation": "int",
            "default": "0",
            "note": "The end site of focus region"
        },
    ],
    "add_heatmap": [
        ["path", "group", "category", "size_factor", "label", "title",
         "barcode_groups", "barcode_tag", "umi_tag", "library", "color",
         "font_size", "show_y_label", "log_trans"],
        {
            "key": "group",
            "annotation": "str",
            "default": "<class 'inspect._empty'>",
            "note": "The heatmap group"
        },
        {
            "key": "do_scale",
            "annotation": "bool",
            "default": "false",
            "note": "Whether to scale the matrix"
        },
        {
            "key": "clustering",
            "annotation": "bool",
            "default": "false",
            "note": "Whether reorder matrix by clustering"
        },
        {
            "key": "clustering_method",
            "annotation": "str",
            "default": "ward",
            "note": "same as  scipy.cluster.hierarchy.linkage"
        },
        {
            "key": "distance_metric",
            "annotation": "str",
            "default": "euclidean",
            "note": "same as scipy.spatial.distance.pdist"
        },
        {
            "key": "show_row_names",
            "annotation": "bool",
            "default": "false",
            "note": "Whether to show row names"
        },
        {
            "key": "vmin",
            "annotation": "float",
            "default": "<class 'inspect._empty'>",
            "note": "Values to anchor the colormap, otherwise they are inferred from the data and other keyword arguments."
        },
        {
            "key": "vmax",
            "annotation": "float",
            "default": "<class 'inspect._empty'>",
            "note": "Values to anchor the colormap, otherwise they are inferred from the data and other keyword arguments."
        }
    ],
    "add_hic": [
        ["path", "category", "label", "color", "font_size", "n_y_ticks", "show_y_label"],
        {
            "key": "trans",
            "annotation": "Optional[str]",
            "default": "<class 'inspect._empty'>",
            "note": "The data transform for HiC matrix, eg log1p, log2, log10 (optional)"
        },
        {
            "key": "show_legend",
            "annotation": "bool",
            "default": "true",
            "note": "Whether to show legend"
        },
        {
            "key": "depth",
            "annotation": "int",
            "default": "30000",
            "note": "The default depth for HiCMatrix"
        },
    ],
    "add_igv": [
        ["path", "category", "label", "library", "font_size", "n_y_ticks", "show_y_label"],
        {
            "key": "exon_focus",
            "annotation": "Optional[str]",
            "default": "<class 'inspect._empty'>",
            "note": "The exon needs to highlight"
        },
        {
            "key": "deletion_ignore",
            "annotation": "Optional[int]",
            "default": "true",
            "note": "Whether ignore the gap in full length sequencing data"
        },
        {
            "key": "del_ratio_ignore",
            "annotation": "float",
            "default": "0.5",
            "note": "Ignore the deletion gap in nanopore or pacbio reads"
        },
        {
            "key": "exon_color",
            "annotation": "color",
            "default": "#000000",
            "note": "The color of exons"
        },
        {
            "key": "intron_color",
            "annotation": "color",
            "default": "#000000",
            "note": "The color of introns"
        },
        {
            "key": "exon_width",
            "annotation": "float",
            "default": "0.3",
            "note": "The default width of exons"
        }
    ],
    "add_interval": [
        {
            "key": "interval",
            "annotation": "str",
            "default": "<class 'inspect._empty'>",
            "note": "Path to interval file"
        },
        {
            "key": "interval_label",
            "annotation": "str",
            "default": "<class 'inspect._empty'>",
            "note": "The label of interval"
        }
    ],
    "add_line": [
        ["path", "group", "category", "size_factor", "label", "title",
         "barcode_groups", "barcode_tag", "umi_tag", "library", "color",
         "font_size", "show_y_label", "log_trans"],
        {
            "key": "show_legend",
            "annotation": "bool",
            "default": "false",
            "note": "Whether to show legend"
        },
        {
            "key": "legend_position",
            "annotation": "str",
            "default": "upper right",
            "note": "The position of legend"
        },
        {
            "key": "legend_ncol",
            "annotation": "int",
            "default": "0",
            "note": "The number of columns of legend"
        }
    ],
    "add_links": [
        {
            "key": "start",
            "annotation": "int",
            "default": "0",
            "note": "The start site of link"
        },
        {
            "key": "end",
            "annotation": "int",
            "default": "0",
            "note": "The end site of link"
        },
        {
            "key": "label",
            "annotation": "str",
            "default": "<class 'inspect._empty'>",
            "note": "The label of link"
        },
        {
            "key": "color",
            "annotation": "color",
            "default": "#409EFF",
            "note": "The color of link"
        }
    ],
    "add_sites": [
        {
            "key": "sites",
            "annotation": "str",
            "default": "",
            "note": "Where to plot additional indicator lines, comma separated int"
        }
    ],
    "add_stroke": [
        {
            "key": "start",
            "annotation": "int",
            "default": "0",
            "note": "The start site of stroke"
        },
        {
            "key": "end",
            "annotation": "int",
            "default": "0",
            "note": "The end site of stroke"
        },
        {
            "key": "label",
            "annotation": "str",
            "default": "<class 'inspect._empty'>",
            "note": "The label of stroke"
        },
        {
            "key": "color",
            "annotation": "color",
            "default": "#000000",
            "note": "The color of stroke"
        }
    ],
    "plot": [
        {
            "key": "reference_scale",
            "annotation": "Union[int, float]",
            "default": "0.25",
            "note": "The size of reference plot in final plot"
        },
        {
            "key": "stroke_scale",
            "annotation": "Union[int, float]",
            "default": "0.25",
            "note": "The size of stroke plot in final image"
        },
        {
            "key": "dpi",
            "annotation": "int",
            "default": "300",
            "note": "The resolution of output file"
        },
        {
            "key": "width",
            "annotation": "Union[int, float]",
            "default": "10",
            "note": "The width of output file, default adjust image width by content"
        },
        {
            "key": "height",
            "annotation": "Union[int, float]",
            "default": "1",
            "note": "The height of single subplot, default adjust image height by content"
        },
        {
            "key": "raster",
            "annotation": "bool",
            "default": "false",
            "note": "The would convert heatmap and site plot to raster image (speed up rendering and produce smaller files), only affects pdf, svg and PS"
        },
        {
            "key": "distance_ratio",
            "annotation": "float",
            "default": "0",
            "note": "The distance between transcript label and transcript line"
        },
        {
            "key": "n_jobs",
            "annotation": "int",
            "default": "1",
            "note": "How many cpu to use"
        },
        {
            "key": "fill_step",
            "annotation": "str",
            "default": "post",
            "note": "Define step if the filling should be a step function, i.e. constant in between x. The value determines where the step will occur:\n"
                    "pre:The y value is continued constantly to the left from every x position, i.e. the interval (x[i-1], x[i]] has the value y[i].\n"
                    "post:The y value is continued constantly to the right from every x position, i.e. the interval [x[i], x[i+1]) has the value y[i].\n"
                    "mid:Steps occur half-way between the x positions."
        },
        {
            "key": "same_y",
            "annotation": "bool",
            "default": "false",
            "note": "Whether different sashimi/line plots shared same y-axis boundaries"
        },
        {
            "key": "remove_duplicate_umi",
            "annotation": "bool",
            "default": "false",
            "note": "Drop duplicated UMIs by barcode"
        },
        {
            "key": "threshold",
            "annotation": "int",
            "default": "0",
            "note": "Threshold to filter low abundance junctions"
        },
        {
            "key": "normalize_format",
            "annotation": "choice['normal', 'cpm', 'rpkm']",
            "default": "normal",
            "note": "The normalize format for bam file"
        },
        {
            "key": "fill_step",
            "annotation": "choice['post', 'pre', 'mid']",
            "default": "post",
            "note": "Define step if the filling should be a step function, i.e. constant in between x. Detailed info please check: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.fill_between.html"
        },
        {
            "key": "normalize_format",
            "annotation": "choice['count', 'cpm', 'rpkm']",
            "default": "count",
            "note": "The normalize format for bam file"
        },
        {
            "key": "smooth_bin",
            "annotation": "int",
            "default": "20",
            "note": "The bin size used to smooth ATAC fragments."
        }
    ],
    "set_reference": [
        {
            "key": "gtf",
            "annotation": "str",
            "default": "<class 'inspect._empty'>",
            "note": "Path to reference file, in sorted and gzipped gtf format"
        },
        {
            "key": "add_domain",
            "annotation": "bool",
            "default": "false",
            "note": "Add domain information into reference track"
        },
        {
            "key": "remove_empty_transcripts",
            "annotation": "bool",
            "default": "false",
            "note": "Whether show transcripts without any exons in target region"
        },
        {
            "key": "choose_primary",
            "annotation": "bool",
            "default": "false",
            "note": "Whether choose primary transcript to plot."
        },
        {
            "key": "color",
            "annotation": "color",
            "default": "#000000",
            "note": "The color of exons"
        },
        {
            "key": "font_size",
            "annotation": "int",
            "default": "5",
            "note": "The font size of reference"
        },
        {
            "key": "no_gene",
            "annotation": "bool",
            "default": "false",
            "note": "Do not show gene id next to transcript id"
        },
        {
            "key": "show_id",
            "annotation": "bool",
            "default": "false",
            "note": "Show gene name or gene id"
        },
        {
            "key": "show_exon_id",
            "annotation": "bool",
            "default": "false",
            "note": "Whether show gene id or gene name"
        },
        {
            "key": "transcripts_to_show",
            "annotation": "str",
            "default": "",
            "note": "Which transcript to show, transcript name or id in gtf file, eg: transcript1,transcript2"
        },
        {
            "key": "intron_scale",
            "annotation": "float",
            "default": "0.5",
            "note": "The scale of intron"
        },
        {
            "key": "exon_scale",
            "annotation": "float",
            "default": "1",
            "note": "The scale of exon"
        },
    ],
    "set_region": [
        {
            "key": "chromosome",
            "annotation": "str",
            "default": "empty",
            "note": "The chromosome of target region"
        },
        {
            "key": "start",
            "annotation": "int",
            "default": "0",
            "note": "The start site of target region"
        },
        {
            "key": "end",
            "annotation": "int",
            "default": "0",
            "note": "The end site of target region"
        },
        {
            "key": "strand",
            "annotation": "str",
            "default": "+",
            "note": "The strand of target region"
        }
    ],
    "set_sequence": [
        {
            "key": "fasta",
            "annotation": "str",
            "default": "<class 'inspect._empty'>",
            "note": "Path to indexed fasta file"
        }
    ]
}


# backend and frontend data format restriction
class Path:
    def __init__(self, path: str, isdir: bool):
        self.path = path
        self.isdir = isdir


class Param:
    def __init__(self, key: str, annotation: str, default: str, note: Optional[str] = None, **kwargs):
        self.key = key
        self.annotation = annotation
        self.default = default
        self.note = note


class PlotParam:

    __slots__ = ["path", "type", "param"]

    def __init__(self, path: str, type: Optional[str], param: List[Param]):
        self.path = path
        self.type = type
        self.param = param

    def __dict__(self):
        return {
            "path": self.path,
            "type": self.type,
            "param": [vars(x) for x in self.param]
        }

    @classmethod
    def create(cls, data: dict):
        if not data:
            return None
        res = {}
        for key, value in data.items():
            if key not in cls.__slots__:
                continue
            if key != "param":
                res[key] = value
            else:
                res[key] = [Param(**val) for val in value]
        return cls(**res)


class PostForm:
    __slots__ = ["region", "reference", "files", "draw"]

    def __init__(self, region: PlotParam, reference: PlotParam, files: List[PlotParam], draw: PlotParam):
        self.region = region
        self.reference = reference
        self.files = files
        self.draw = draw

    @classmethod
    def create(cls, data: dict):
        res = {}
        for key, value in data.items():
            if key not in cls.__slots__:
                continue
            if key != "files":
                if value:
                    res[key] = PlotParam.create(value)
            else:
                res[key] = [PlotParam.create(val) for val in value if val]

        return cls(**res)


class Logs:
    def __init__(self, time: str, level: str, source: str, message: str):
        self.time = time
        self.level = level
        self.source = source
        self.message = message


def get_value(param: Param):
    u"""
    convert parameters in post param to specific format
    """

    if "empty" in str(param.default):
        return None

    if "str" in param.annotation:
        return param.default if str(param.default).lower() != "false" else None

    if param.annotation == "int":
        return int(param.default)

    if param.annotation == "float":
        return float(param.default)

    if param.annotation == "bool":
        return str(param.default).lower() == "true"

    if param.annotation == "Union[int, float]":
        return float(param.default)

    if param.annotation.startswith("choice"):
        return param.default.lower()

    return param.default


def clean(annotation):
    if "empty" in annotation:
        return ""

    annotation = annotation.replace("<class '", "").replace("'>", "").replace("'", "")
    return annotation.replace("typing.", "")


@app.route("/")
def root():
    return render_template("index.html")


@app.route("/api/file")
def file():
    target = request.args.get('target', __DIR__)
    valid = request.args.get("valid", "false").lower() == "true"
    fs = []
    src = target

    if valid:
        return jsonify(os.path.isfile(target))

    try:
        if os.path.isdir(src):
            src = os.path.join(src, "*")
        else:
            src += "*"
        for x in glob(src):
            if not os.path.basename(x).startswith("."):
                fs.append(Path(path=x, isdir=os.path.isdir(x)))
    except OSError as err:
        # raise HTTPException(status_code=404, detail=str(err))
        return jsonify(str(err)), 404

    return jsonify([vars(x) for x in sorted(fs, key=lambda x: (not x.isdir, x.path))])


@app.route("/api/params")
def params():
    target = request.args.get("target")
    res = []

    for p in __PARAMS__[target]:
        if not isinstance(p, dict):
            for cp in __COMMON_PARAMS__:
                if cp["key"] in p:
                    res.append(cp)
        else:
            res.append(p)

    return jsonify(res)


@app.route("/api/plot/<pid>", methods=["POST"])
def plot(pid: str):
    param = request.get_json()
    param = PostForm.create(param)

    log = os.path.join(__PLOT__, pid + ".log")

    if os.path.exists(log):
        os.remove(log)

    p = Plot(log)

    if param:
        with open(os.path.join(__PLOT__, pid), "wb+") as w:
            pickle.dump(param, w)
    else:
        with open(os.path.join(__PLOT__, pid), "rb") as r:
            param = pickle.load(r)

    def plot_form_to_dict(param_: PlotParam):
        kwargs = {}
        for para_ in param_.param:
            val = get_value(para_)
            if val is not None:
                kwargs[para_.key] = val

        if "no_gene" in kwargs:
            kwargs["show_gene"] = not kwargs.pop("no_gene")

        if "distance_ratio" in kwargs:
            kwargs["distance_between_label_axis"] = kwargs["distance_ratio"]
        return kwargs

    try:
        # set region
        p.set_region(**plot_form_to_dict(param.region))

        # set reference
        p.set_reference(param.reference.path, **plot_form_to_dict(param.reference))

        # set files
        for param_ in param.files:
            kwargs = plot_form_to_dict(param_)
            if re.search(r"_(sites|stroke|link|focus)", param_.type.lower()):
                getattr(p, param_.type.lower())(**kwargs)
            else:
                if "barcode_groups" in kwargs.keys() and kwargs["barcode_groups"]:
                    barcodes, sc_colors = load_barcodes(kwargs["barcode_groups"])
                    kwargs["barcode_groups"] = barcodes
                    for group in barcodes[kwargs["label"]].keys():
                        kwargs["barcode"] = group
                        kwargs["color"] = sc_colors[group]
                        kwargs["y_label"] = group
                        getattr(p, param_.type.lower())(param_.path, **kwargs)
                else:
                    for tag in ["vmin", "vmax"]:
                        if tag in kwargs.keys():
                            try:
                                kwargs[tag] = float(kwargs[tag])
                            except Exception:
                                kwargs[tag] = None

                    getattr(p, param_.type.lower())(param_.path, **kwargs)

        # draw
        param_ = param.draw
        kwargs = plot_form_to_dict(param_)
        if param_.type == "plot":
            o = p.plot(return_image="png", **kwargs)
            o.seek(0)
            return send_file(o, mimetype=f"image/png")
        elif param_.type == "save":
            o = p.plot(return_image="pdf", **kwargs)
            o.seek(0)
            return send_file(o, mimetype=f"image/pdf", as_attachment=True, download_name=f"{p.region}.pdf")

    except Exception as err:
        logger.error(err)
        return jsonify(str(err)), 501
    return jsonify("")

@app.route("/api/del")
def delete():
    pid = request.args.get("pid", "?")
    pk = os.path.join(__PLOT__, pid)
    if os.path.exists(pk):
        os.remove(pk)

    log = os.path.join(__PLOT__, pid + ".log")
    if os.path.exists(log):
        os.remove(log)
    return jsonify("done")


@app.route("/api/log")
def logs():
    pid = request.args.get("pid", "?")
    debug = request.args.get("debug", "false").lower() == "true"
    download = request.args.get("download", "false").lower() == "true"

    logfile = os.path.join(__PLOT__, pid + ".log")

    if download:
        return send_file(logfile, as_attachment=True, download_name=os.path.basename(logfile))

    log_info = []
    if os.path.exists(logfile):
        with open(logfile) as r:
            for line in r:
                try:
                    time, level, info = line.strip().split("|")
                    infos = info.split("-")
                    source = infos[0]

                    if not debug and "DEBUG" in level:
                        continue

                    log_info.append(vars(Logs(
                        time=time.strip(),
                        level=level.strip(),
                        source=source.strip(),
                        message=info.replace(source, "").strip()
                    )))
                except Exception as err:
                    log_info.append(vars(Logs(
                        time="",
                        level="WARN",
                        source="log api",
                        message=str(err)
                    )))

    return jsonify(log_info[::-1])


@click.command(context_settings=dict(help_option_names=['--help']), )
@click.option("-h", "--host", type=click.STRING, default="127.0.0.1", help="the ip address binding to")
@click.option("-p", "--port", type=click.INT, default=5000, help="the port to listen on")
@click.option("--plots", type=click.Path(), default=__PLOT__,
              help="the path to directory where to save the backend plot data and logs, required while using appImage.")
@click.version_option(__version__, message="Current version %(version)s")
def main(host: str, port: int, plots: str):
    global __PLOT__
    if plots:
        __PLOT__ = plots
    os.makedirs(__PLOT__, exist_ok=True)

    app.run(host=host, port=port)


if __name__ == '__main__':
    main()

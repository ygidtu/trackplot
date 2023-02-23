#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
Web UI of sashimi
"""

import pickle
import sys
import traceback
from glob import glob

import click
import uvicorn
from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel

from sashimi.cli import load_barcodes
from sashimi.plot import *

__DIR__ = os.path.abspath(os.path.dirname(__file__))
__UI__ = os.path.join(__DIR__, "ui")

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.mount("/static", StaticFiles(directory=__UI__, html=True), name="static")
templates = Jinja2Templates(directory=__UI__)

__PLOT__ = os.path.join(os.path.dirname(__file__), "plots")
os.makedirs(__PLOT__, exist_ok=True)

__COMMON_PARAMS__ = [
    {
        "key": "path",
        "annotation": "str",
        "default": "<class 'inspect._empty'>",
        "note": "Please input the path to input file"
    },
    {
        "key": "category",
        "annotation": "str",
        "default": "bam",
        "note": "The category of input input file"
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
        "key": "barcode_tag",
        "annotation": "str",
        "default": "BC",
        "note": "The barcode tag in bam file"
    },
    {
        "key": "umi_tag",
        "annotation": "str",
        "default": "UB",
        "note": "The umi tag in bam file"
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
        "default": "True",
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
        "annotation": "",
        "default": "blue",
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
            "default": "False",
            "note": "Whether draw density plot by strand"
        },
        {
            "key": "show_junction_number",
            "annotation": "bool",
            "default": "True",
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
            "default": "False",
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
            "default": "False",
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
            "default": "False",
            "note": "Whether to scale the matrix"
        },
        {
            "key": "clustering",
            "annotation": "bool",
            "default": "False",
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
            "default": "False",
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
            "default": "True",
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
            "default": "True",
            "note": ""
        },
        {
            "key": "del_ratio_ignore",
            "annotation": "float",
            "default": "0.5",
            "note": ""
        },
        {
            "key": "exon_color",
            "annotation": "Optional[str]",
            "default": "black",
            "note": "The color of exons"
        },
        {
            "key": "intron_color",
            "annotation": "Optional[str]",
            "default": "black",
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
            "default": "False",
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
            "annotation": "str",
            "default": "blue",
            "note": "The color of link"
        }
    ],
    "add_sites": [
        {
            "key": "sites",
            "annotation": "",
            "default": "<class 'inspect._empty'>",
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
            "annotation": "str",
            "default": "black",
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
            "default": "0",
            "note": "The width of output file, default adjust image width by content"
        },
        {
            "key": "height",
            "annotation": "Union[int, float]",
            "default": "0",
            "note": "The height of single subplot, default adjust image height by content"
        },
        {
            "key": "raster",
            "annotation": "bool",
            "default": "False",
            "note": "The would convert heatmap and site plot to raster image (speed up rendering and produce smaller files), only affects pdf, svg and PS"
        },
        {
            "key": "distance_between_label_axis",
            "annotation": "float",
            "default": "0.3",
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
            "note": "Define step if the filling should be a step function, i.e. constant in between x. "
                    "Detailed info please check: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.fill_between.html"
        },
        {
            "key": "same_y",
            "annotation": "bool",
            "default": "False",
            "note": "Whether different sashimi/line plots shared same y-axis boundaries"
        },
        {
            "key": "remove_duplicate_umi",
            "annotation": "bool",
            "default": "False",
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
            "default": "False",
            "note": "Add domain information into reference track"
        },
        {
            "key": "local_domain",
            "annotation": "Optional[str]",
            "default": "False",
            "note": "Load local domain folder and load into reference track, download from "
                    "https://hgdownload.soe.ucsc.edu/gbdb/hg38/uniprot/"
        },
        {
            "key": "interval",
            "annotation": "Optional[str]",
            "default": "<class 'inspect._empty'>",
            "note": "Path to list of interval files in bed format, 1st column is path to file, 2nd column is the label [optional]"
        },
        {
            "key": "interval_label",
            "annotation": "Optional[str]",
            "default": "<class 'inspect._empty'>",
            "note": "The label of interval"
        },
        {
            "key": "transcripts",
            "annotation": "Optional[List[str]]",
            "default": "<class 'inspect._empty'>",
            "note": "The name transcripts for "
        },
        {
            "key": "remove_empty_transcripts",
            "annotation": "bool",
            "default": "False",
            "note": "Whether show transcripts without any exons in target region"
        },
        {
            "key": "choose_primary",
            "annotation": "bool",
            "default": "False",
            "note": "Whether choose primary transcript to plot."
        },
        {
            "key": "color",
            "annotation": "Optional[str]",
            "default": "black",
            "note": "The color pf transcript"
        },
        {
            "key": "font_size",
            "annotation": "int",
            "default": "5",
            "note": "The font size of reference"
        },
        {
            "key": "show_gene",
            "annotation": "bool",
            "default": "False",
            "note": "Show gene name or gene id"
        },
        {
            "key": "show_id",
            "annotation": "bool",
            "default": "False",
            "note": "Show gene id/transcript id instead of gene name/transcript name"
        },
        {
            "key": "exon_width",
            "annotation": "float",
            "default": "0.3",
            "note": "The width of exons"
        },
        {
            "key": "show_exon_id",
            "annotation": "bool",
            "default": "False",
            "note": "Show exon id"
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


class Path(BaseModel):
    path: str
    isdir: bool


class Param(BaseModel):
    key: str
    annotation: str
    default: str
    note: Optional[str]


class PlotParam(BaseModel):
    path: str
    param: List[Param]


def get_value(param: Param):
    if "empty" in param.default:
        return None
    if "str" in param.annotation:
        return param.default if param.default != "False" else None

    if param.annotation == "int":
        return int(param.default)

    if param.annotation == "float":
        return float(param.default)

    if param.annotation == "bool":
        return param.default.lower() != 'false'

    if param.annotation == "Union[int, float]":
        return float(param.default)

    if param.annotation.startswith("choice"):
        return param.default.lower()

    return param.default


@app.get("/")
async def root(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.get("/api/file")
async def file(target: str = __DIR__, valid: bool = False) -> Union[bool, List[Path]]:
    fs = []
    src = target

    if valid and target:
        return os.path.isfile(target)

    try:
        if os.path.isdir(src):
            src = os.path.join(src, "*")
        else:
            src += "*"
        for x in glob(src):
            if not os.path.basename(x).startswith("."):
                fs.append(Path(path=x, isdir=os.path.isdir(x)))
    except OSError as err:
        raise HTTPException(status_code=404, detail=str(err))

    return sorted(fs, key=lambda x: (not x.isdir, x.path))


def clean(annotation):
    if "empty" in annotation:
        return ""

    annotation = annotation.replace("<class '", "").replace("'>", "").replace("'", "")
    return annotation.replace("typing.", "")


@app.get("/api/params")
async def params(target: Optional[str] = None):
    res = []

    for p in __PARAMS__[target]:
        if not isinstance(p, dict):
            for cp in __COMMON_PARAMS__:
                if cp["key"] in p:
                    res.append(Param(**cp))
        else:
            res.append(Param(**p))

    # if not target:
    #     method_list = [attribute for attribute in dir(Plot) if
    #                    callable(getattr(Plot, attribute)) and attribute.startswith('__') is False]
    #     return method_list
    #
    # sig = inspect.signature(getattr(Plot, target))
    #
    # for param in sig.parameters.values():
    #     if param.name == "self":
    #         continue
    #
    #     if param.name == "output":
    #         continue
    #
    #     if param.name == "return_image":
    #         continue
    #
    #     if param.name == "sc_height_ratio":
    #         continue
    #
    #     if param.name == "barcode":
    #         continue
    #
    #     res.append(Param(
    #         key=param.name,
    #         annotation=clean(str(param.annotation)),
    #         default=str(param.default) if param.default is not None else ""
    #     ))
    #
    # if target == "plot":
    #     res.append(Param(key="same_y", annotation='bool', default="False"))
    #     res.append(Param(key="remove_duplicate_umi", annotation='bool', default="False"))
    #     res.append(Param(key="threshold", annotation='int', default="0"))
    #     res.append(Param(key="normalize_format", annotation="choice['normal', 'cpm', 'rpkm']",
    #                      default="normal", note="The normalize format for bam file"))
    #     res.append(Param(key="fill_step", annotation="choice['post', 'pre', 'mid']", default="post",
    #                      note="Define step if the filling should be a step function, i.e. constant in between x. "
    #                           "Detailed info please check: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.fill_between.html"))
    #     res.append(Param(key="smooth_bin", annotation='int',
    #                      default="20", note="The bin size used to smooth ATAC fragments."))

    return res


@app.post("/api/plot")
async def plot(pid: str, param: PlotParam, func: str):
    pk = os.path.join(__PLOT__, pid)
    if not os.path.exists(pk):
        p = Plot()
        with open(pk, "wb+") as w:
            pickle.dump(p, w)
    else:
        with open(pk, "rb") as r:
            p = pickle.load(r)

    kwargs = {}
    for param_ in param.param:
        val = get_value(param_)
        if val is not None:
            kwargs[param_.key] = val

    try:
        if func == "plot":
            o = p.plot(return_image="png", **kwargs)
            o.seek(0)
            return StreamingResponse(o, media_type=f"image/png")
        elif func == "save":
            o = p.plot(return_image="pdf", **kwargs)
            o.seek(0)
            resp = StreamingResponse(o, media_type=f"image/pdf")
            resp.headers["Access-Control-Expose-Headers"] = "Content-Disposition"
            resp.headers["Content-Disposition"] = str(p.region) + ".pdf"
            return resp
        elif func == "set_region":
            p.set_region(**kwargs)
        else:
            if "barcode_groups" in kwargs.keys() and kwargs["barcode_groups"]:
                barcodes, sc_colors = load_barcodes(kwargs["barcode_groups"])
                kwargs["barcode_groups"] = barcodes
                for group in barcodes[kwargs["label"]].keys():
                    kwargs["barcode"] = group
                    kwargs["color"] = sc_colors[group]
                    kwargs["y_label"] = group
                    getattr(p, func)(param.path, **kwargs)
            else:
                for tag in ["vmin", "vmax"]:
                    if tag in kwargs.keys():
                        try:
                            kwargs[tag] = float(kwargs[tag])
                        except Exception:
                            kwargs[tag] = None

                getattr(p, func)(param.path, **kwargs)
    except (AssertionError, OSError, TypeError) as err:
        print(err)
        exc_type, exc_value, exc_traceback = sys.exc_info()
        raise HTTPException(status_code=404,
                            detail=f"{err}: {traceback.format_exception(exc_type, exc_value, exc_traceback)}")

    with open(pk, "wb+") as w:
        pickle.dump(p, w)


@app.get("/api/del")
async def delete(pid: str):
    pk = os.path.join(__PLOT__, pid)
    if os.path.exists(pk):
        os.remove(pk)


@click.command(context_settings=dict(help_option_names=['--help']), )
@click.option("-h", "--host", type=click.STRING, default="127.0.0.1", help="the ip address binding to")
@click.option("-p", "--port", type=click.INT, default=5000, help="the port to listen on")
@click.option("-r", "--reload", is_flag=True, help="auto-reload for development")
def main(host: str, port: int, reload: bool):
    uvicorn.run(app='server:app', host=host, port=int(port), log_level="info", reload=reload)


if __name__ == '__main__':
    main()

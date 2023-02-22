#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
Web UI of sashimi
"""

import inspect
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


class Path(BaseModel):
    path: str
    isdir: bool


class Param(BaseModel):
    key: str
    annotation: str
    default: str


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

    if not target:
        method_list = [attribute for attribute in dir(Plot) if
                       callable(getattr(Plot, attribute)) and attribute.startswith('__') is False]
        return method_list

    sig = inspect.signature(getattr(Plot, target))

    for param in sig.parameters.values():
        if param.name == "self":
            continue

        if param.name == "output":
            continue

        if param.name == "return_image":
            continue

        if param.name == "sc_height_ratio":
            continue

        if param.name == "barcode":
            continue

        res.append(Param(
            key=param.name,
            annotation=clean(str(param.annotation)),
            default=str(param.default) if param.default is not None else ""
        ))

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
def main(host: str, port: int):
    uvicorn.run(app='server:app', host=host, port=int(port), log_level="info", reload=False)


if __name__ == '__main__':
    main()

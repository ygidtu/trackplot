#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
Web UI of sashimi
"""
import pickle
import re
from glob import glob

from flask import Flask, render_template, jsonify, send_file, request

from trackplot.conf.ui import __SUPPORT_FORMAT__, __COMMON_PARAMS__, __PARAMS__
from trackplot.plot import *


__DIR__ = os.path.abspath(os.path.dirname(__file__))
__UI__ = os.path.join(__DIR__, "../ui")
__PLOT__ = os.path.join(os.path.dirname(__file__), "../plots")


app = Flask(__name__, static_url_path="/static", static_folder=__UI__, template_folder=__UI__)

# from flask_cors import CORS
# CORS(app)

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
    __slots__ = ["region", "annotation", "files", "draw"]

    def __init__(self, region: PlotParam, annotation: PlotParam, files: List[PlotParam], draw: PlotParam):
        self.region = region
        self.annotation = annotation
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
        return float(param.default) if param.default else 0.0

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
    target = request.args.get('target', app.config["DATA"])

    if not target.strip("/"):
        target = app.config["DATA"]

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

    category = None
    for p in __PARAMS__[target]:
        if not isinstance(p, dict):
            for cp in __COMMON_PARAMS__:
                if cp["key"] in p:

                    if cp["key"] == "category" and target in __SUPPORT_FORMAT__:
                        category = cp
                        category["annotation"] = f"choice[{','.join(__SUPPORT_FORMAT__[target])}]"
                        category["default"] = __SUPPORT_FORMAT__[target][0]
                    else:
                        res.append(cp)
        else:
            res.append(p)

    if category:
        res.append(category)

    return jsonify(res)


@app.route("/api/plot/<pid>", methods=["POST"])
def plot(pid: str):
    param = request.get_json()
    param = PostForm.create(param)

    log = os.path.join(app.config["PLOTTING"], pid + ".log")

    if os.path.exists(log):
        os.remove(log)

    p = Plot(logfile=log, backend="agg")

    if param:
        with open(os.path.join(app.config["PLOTTING"], pid), "wb+") as w:
            pickle.dump(param, w)
    else:
        with open(os.path.join(app.config["PLOTTING"], pid), "rb") as r:
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

        # set annotation
        p.set_annotation(param.annotation.path, **plot_form_to_dict(param.annotation))

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
        logger.exception(err)
        return jsonify(str(err)), 501
    return jsonify("")


@app.route("/api/del")
def delete():
    pid = request.args.get("pid", "?")
    pk = os.path.join(app.config["PLOTTING"], pid)
    if os.path.exists(pk):
        os.remove(pk)

    log = os.path.join(app.config["PLOTTING"], pid + ".log")
    if os.path.exists(log):
        os.remove(log)
    return jsonify("done")


@app.route("/api/log")
def logs():
    pid = request.args.get("pid", "?")
    debug = request.args.get("debug", "false").lower() == "true"
    download = request.args.get("download", "false").lower() == "true"

    logfile = os.path.join(app.config["PLOTTING"], pid + ".log")

    if download:
        return send_file(logfile, as_attachment=True, download_name=os.path.basename(logfile))

    log_info = []
    if os.path.exists(logfile):
        with open(logfile) as r:
            for line in r:
                try:
                    if "|" not in line and len(log_info) > 0:
                        log_info[-1]['message'] = f"{log_info[-1]['message']}\n{line}"
                    else:
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


def run(host: str, port: int, plots: str, data: str):
    os.makedirs(plots, exist_ok=True)
    app.config["PLOTTING"] = plots
    app.config["DATA"] = data

    app.run(host=host, port=port, debug=False)


if __name__ == '__main__':
    pass

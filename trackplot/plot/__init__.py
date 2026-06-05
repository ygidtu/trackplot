#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
trackplot.plot — plot orchestration and rendering package.

Re-exports all public symbols for backward compatibility.
Use ``from trackplot.plot.core import Plot`` or ``from trackplot.plot.render import plot_density`` for direct access.
"""

from trackplot.plot.coord import (init_graph_coords, set_focus,
                                  set_indicator_lines, set_x_ticks,
                                  set_y_ticks)
from trackplot.plot.core import (Plot, __author__, __email__, __version__,
                                 add_object_error)
from trackplot.plot.info import PlotInfo
from trackplot.plot.limits import precompute_y_limits
from trackplot.plot.render import (plot_annotation, plot_density, plot_heatmap,
                                   plot_hic, plot_igv_like, plot_line,
                                   plot_links, plot_motif, plot_site_plot,
                                   plot_stroke)
from trackplot.plot.utils import (__merge_exons__, cubic_bezier,
                                  get_limited_index, load_barcodes,
                                  make_text_elements)

__all__ = [
    "Plot",
    "PlotInfo",
    "add_object_error",
    "__version__",
    "__author__",
    "__email__",
    # utils
    "load_barcodes",
    "get_limited_index",
    "cubic_bezier",
    "__merge_exons__",
    "make_text_elements",
    # coord
    "init_graph_coords",
    "set_x_ticks",
    "set_y_ticks",
    "set_focus",
    "set_indicator_lines",
    # limits
    "precompute_y_limits",
    # render
    "plot_stroke",
    "plot_annotation",
    "plot_density",
    "plot_site_plot",
    "plot_heatmap",
    "plot_hic",
    "plot_line",
    "plot_igv_like",
    "plot_links",
    "plot_motif",
]

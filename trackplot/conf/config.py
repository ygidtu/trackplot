#!/usr/bin/env python3
# -*- coding:utf-8 -*-

IMAGE_TYPE = ["annotation", "interval", "density", "heatmap", "line", "igv", "hic", "motif"]

COLORS = ["#A0C807", "#0084D1", '#db5f57', '#dbc257', '#91db57', '#57db80', '#57d3db', '#5770db', '#a157db', '#db57b2']
COLORMAP = [
    'Blues', 'Greens', 'Oranges', 'Reds', 'Greys', 'Purples',
    'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
    'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
    'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
    'viridis', 'plasma', 'inferno', 'magma', 'cividis',
    'pink', 'spring', 'summer', 'autumn', 'winter',
    'cool', 'Wistia', 'hot', 'afmhot', 'gist_heat',
    'copper', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
    'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic', "RdYlBu_r"
]

DISTANCE_METRIC = [
    "braycurtis", "canberra", "chebyshev",
    "cityblock", "correlation", "cosine",
    "dice", "euclidean", "hamming",
    "jaccard", "jensenshannon", "kulsinski",
    "kulczynski1", "mahalanobis", "matching",
    "minkowski", "rogerstanimoto", "russellrao",
    "seuclidean", "sokalmichener", "sokalsneath",
    "sqeuclidean", "yule"
]

CLUSTERING_METHOD = ["single", "complete", "average", "weighted", "centroid", "median", "ward"]

NORMALIZATION = ["count", "cpm", "rpkm"]


if __name__ == '__main__':
    pass

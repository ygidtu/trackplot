#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Web UI config
"""
from trackplot.conf.config import COLORMAP


__SUPPORT_FORMAT__ = {
    "add_density": ['bam', 'bigwig', 'bedgraph', 'depth', 'atac'],
    "add_line": ['bam', 'bigwig', 'bedgraph', 'depth', 'atac'],
    "add_heatmap": ['bam', 'bigwig', 'bedgraph', 'depth', 'atac'],
    "add_igv": ['bam', 'bed'],
    "add_hic": ['h5'],

}

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
        "annotation": "choice",
        "default": "bam",
        "note": "The input file category can be selected from options such as BAM, ATAC, IGV, Hi-C, BigWig/BW, BedGraph/BG, and Depth."
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
        "note": "The strand pannotation of input file, frf => fr-firststrand; frs => fr-secondstrand; fru => fr-unstrand"
    },
    {
        "key": "n_y_ticks",
        "annotation": "int",
        "default": "4",
        "note": "Set the number of ticks for the y-axis."
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
        "note": "The color of the current track."
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
        "note": "Choose whether to perform a log transformation for coverage. Use '0' for no transformation, '2' for log2 transformation, and '10' for log10 transformation."
    },
    {
        "key": "group",
        "annotation": "str",
        "default": "default",
        "note": "The group name of heatmap/line track"
    },
]


__PARAMS__ = {
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
            "note": "Whether draw density plot in a strand-aware manner. If enabled, the tool will consider the library parameter (frf or frs) to perform strand-aware coverage calculation."
        },
        {
            "key": "show_junction_number",
            "annotation": "bool",
            "default": "true",
            "note": "Whether to show the depth of each splicing junction. If activated, the tool will display the depth of each splicing junction which will be visualized at the midpoint of the junction span line."
        },
        {
            "key": "junction_number_font_size",
            "annotation": "int",
            "default": "5",
            "note": "The font size of the depth of splicing junction."
        },
        {
            "key": "show_site_plot",
            "annotation": "bool",
            "default": "false",
            "note": "Whether to include a site coverage plot for the current track as well. If enabled, only the most 3' end site of each read will be counted in the coverage matrix. When strand-aware mode is enabled, this allows for easy distinction of alternative polyadenylation events in 3'-tag sequencing datasets."
        },
        {
            "key": "strand_choice",
            "annotation": "choice['all', '+', '-']",
            "default": "all",
            "note": "The coverage matrix used for generating the site plot can be specified by strand. By default, coverage from all strands is used. If the + or - option is enabled, only the coverage matrix from the chosen strand will be displayed."
        }
    ],
    "add_focus": [
        {
            "key": "start",
            "annotation": "int",
            "default": "0",
            "note": "The start site of focus region. More detailed instructions please refer to https://trackplot.readthedocs.io/en/latest/command/#additional-annotation"
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
         "barcode_groups", "barcode_tag", "umi_tag", "library",
         "font_size", "show_y_label", "log_trans"],
        {
            "key": "color",
            "annotation": f"select[{','.join(COLORMAP)}]",
            "default": "RdYlBu_r",
            "note": "The color of the current track."
        },
        {
            "key": "do_scale",
            "annotation": "bool",
            "default": "false",
            "note": "Whether perform a scale for the coverage matrix by samples."
        },
        {
            "key": "clustering",
            "annotation": "bool",
            "default": "false",
            "note": "Whether perform a cluster analysis for the samples."
        },
        {
            "key": "clustering_method",
            "annotation": "str",
            "default": "ward",
            "note": "If `clustering` enabled, user could select a clustering method for downstream analysis. More clustering method please refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html ."
        },
        {
            "key": "distance_metric",
            "annotation": "str",
            "default": "euclidean",
            "note": "The method to perform pairwise distances analysis. More detail please refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html ."
        },
        {
            "key": "show_row_names",
            "annotation": "bool",
            "default": "false",
            "note": "Whether to add the filename as the row name of the heatmap."
        },
        {
            "key": "vmin",
            "annotation": "float",
            "default": "<class 'inspect._empty'>",
            "note": "The parameter Vmin determines the minimum value used for the colormap. If not specified, the minimum value will be inferred from the data and other keyword arguments."
        },
        {
            "key": "vmax",
            "annotation": "float",
            "default": "<class 'inspect._empty'>",
            "note": "The parameter Vmax determines the maximum value used for the colormap. If not specified, the maximum value will be inferred from the data and other keyword arguments."
        }
    ],
    "add_hic": [
        ["path", "category", "label", "font_size", "n_y_ticks", "show_y_label"],
        {
            "key": "color",
            "annotation": f"select[{','.join(COLORMAP)}]",
            "default": "RdYlBu_r",
            "note": "The color of the current track."
        },
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
            "note": "Whether to show the legend of each line."
        },
        {
            "key": "depth",
            "annotation": "int",
            "default": "30000",
            "note": "The length of region to visualize genomic interaction."
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
            "note": "Whether to ignore deletion regions in full-length sequencing data can be specified. Full-length sequencing data, especially nanopore data, often contains a high ratio of deletions. By enabling this option, users can exclude these deletion regions from visualization."
        },
        {
            "key": "del_ratio_ignore",
            "annotation": "float",
            "default": "0.5",
            "note": "The threshold for ignoring deletions in full-length sequencing data. For example, if a deletion ignore ratio of 0.5 is selected and a 1000 bp isoform is captured in sequencing, deletions within 1000 * 0.5 length will be filled and not visualized."
        },
        {
            "key": "exon_color",
            "annotation": "color",
            "default": "#000000",
            "note": "The color of exon structure in current track."
        },
        {
            "key": "intron_color",
            "annotation": "color",
            "default": "#000000",
            "note": "The color of intron structure in current track."
        },
        {
            "key": "exon_width",
            "annotation": "float",
            "default": "0.3",
            "note": "The exon width of current track"
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
            "note": "Whether to show the legend of each line."
        },
        {
            "key": "legend_position",
            "annotation": "str",
            "default": "upper right",
            "note": "The position of the legend for current track, default: upper right. Other layout such as upper left, lower left, lower right, right, center left, center right, lower center, upper center, and center will place the legend at the corresponding corner of the track."
        },
        {
            "key": "legend_ncol",
            "annotation": "int",
            "default": "0",
            "note": "The number of columns for current legend."
        }
    ],
    "add_links": [
        {
            "key": "start",
            "annotation": "int",
            "default": "0",
            "note": "The start site of link line. If the \"Links\" option is enabled, a pseudo span line will be added at the bottom of the tracks. This can be used by users to indicate the back-splice junction of circular RNA (circRNA). A detailed instruction please refer to https://trackplot.readthedocs.io/en/latest/command/#circrna-plot ."
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
            "default": "",
            "note": "The label of link"
        },
        {
            "key": "color",
            "annotation": "color",
            "default": "#409EFF",
            "note": "The color of link line."
        }
    ],
    "add_sites": [
        {
            "key": "sites",
            "annotation": "str",
            "default": "",
            "note": "Where to plot additional indicator lines. The user is required to provide an integer number within the region of interest. If there are multiple lines, the integer numbers should be separated by commas."
        }
    ],
    "add_stroke": [
        {
            "key": "start",
            "annotation": "int",
            "default": "0",
            "note": "The start site of stroke. More detailed instructions please refer to https://trackplot.readthedocs.io/en/latest/command/#additional-annotation"
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
            "key": "annotation_scale",
            "annotation": "Union[int, float]",
            "default": "0.25",
            "note": "The height of annotation track, the height = (number of annotation) * annotation_scale"
        },
        {
            "key": "stroke_scale",
            "annotation": "Union[int, float]",
            "default": "0.25",
            "note": "The height of stroke track, the height = (number of strokes) * stroke_scale"
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
            "note": "The width of output file"
        },
        {
            "key": "height",
            "annotation": "Union[int, float]",
            "default": "1",
            "note": "The height of each track, the igv/annotation/stroke tracks will adjust the height according to their scales"
        },
        {
            "key": "raster",
            "annotation": "bool",
            "default": "false",
            "note": "Choose whether to convert the heatmap and site plot to a raster image format, which can speed up rendering and result in smaller file sizes. Note that this option only affects the output formats of PDF, SVG, and PS. If enabled, the resolution of the output will decrease."
        },
        {
            "key": "distance_ratio",
            "annotation": "float",
            "default": "0",
            "note": "Adjust the distance between transcript label and transcript structure in the annotation track. And the distacne between track title and the track"
        },
        {
            "key": "n_jobs",
            "annotation": "int",
            "default": "1",
            "note": "Number of cores to perform data preparation."
        },
        {
            "key": "fill_step",
            "annotation": "str",
            "default": "post",
            "note": "The step parameter of matplotlib.axes.Axes.fill_between. Define step if the filling should be a step function, i.e. constant in between x. The value determines where the step will occur:\n"
                    "pre:The y value is continued constantly to the left from every x position, i.e. the interval (x[i-1], x[i]] has the value y[i].\n"
                    "post:The y value is continued constantly to the right from every x position, i.e. the interval [x[i], x[i+1]) has the value y[i].\n"
                    "mid:Steps occur half-way between the x positions."
        },
        {
            "key": "same_y",
            "annotation": "bool",
            "default": "false",
            "note": "13.Whether share a same y-axis among all data tracks."
        },
        {
            "key": "threshold",
            "annotation": "int",
            "default": "0",
            "note": "The threshold for filter these splicing junction by depth."
        },
        {
            "key": "normalize_format",
            "annotation": "choice['normal', 'cpm', 'rpkm']",
            "default": "normal",
            "note": "The normalize format for bam file."
        },
        {
            "key": "smooth_bin",
            "annotation": "int",
            "default": "20",
            "note": "The bin size used to smooth ATAC fragments."
        }
    ],
    "set_annotation": [
        {
            "key": "gtf",
            "annotation": "str",
            "default": "<class 'inspect._empty'>",
            "note": "Path to annotation file, in sorted and gzipped gtf format"
        },
        {
            "key": "add_domain",
            "annotation": "bool",
            "default": "false",
            "note": "Whether add domain information into the annotation track at bottom. "
                    "When the parameter is activated, "
                    "Trackplot will initiate a request to the EBI and UNIPROT APIs to retrieve protein annotation. "
                    "The retrieved annotation is then converted into genomic coordinates, "
                    "and the resulting domain information is visualized in a manner similar to a transcript structure. "
                    "This can be seen in Supplementary Figure 2A."
        },
        {
            "key": "remove_empty_transcripts",
            "annotation": "bool",
            "default": "false",
            "note": "Whether show transcripts without any exons in the target region. "
                    "Sometimes, "
                    "the region of interest may be located within a large intron of a transcript, "
                    "which contains less informative data. "
                    "To exclude these transcript structures and focus on the relevant regions, "
                    "users can activate this parameter."
        },
        {
            "key": "choose_primary",
            "annotation": "bool",
            "default": "false",
            "note": "Whether choose primary transcript to plot. "
                    "Once the parameter is activated, "
                    "the tool will display the longest transcript of each gene within the region of interest."
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
            "note": "The font size of the annotation tracks."
        },
        {
            "key": "no_gene",
            "annotation": "bool",
            "default": "false",
            "note": "Do not show gene id next to transcript id. "
                    "If set to false, the format \"isoform1|gene_id\" will be displayed. "
                    "If set to true, only the \"isoform1\" ID will be shown."
        },
        {
            "key": "show_id",
            "annotation": "bool",
            "default": "false",
            "note": "Show gene name (ensemble ID) or gene id (gene symbol). "
                    "If set to false, an ensemble ID will be presented. "
                    "If set to true, a symbol ID will be shown."
        },
        {
            "key": "show_exon_id",
            "annotation": "bool",
            "default": "false",
            "note": "Show exon id/exon name."
        },
        {
            "key": "transcripts_to_show",
            "annotation": "str",
            "default": "",
            "note": "To display specific transcripts, "
                    "their names or IDs must be included in the annotation file (GTF file). "
                    "If there are multiple transcripts to show, "
                    "please separate each transcript by a comma (iso1, iso2, iso3)."
        },
        {
            "key": "intron_scale",
            "annotation": "float",
            "default": "0.5",
            "note": "The degree of intron scaling. "
                    "The intron shrinkage scale determines the degree of intron length reduction. "
                    "A smaller value will result in a shorter intron length,"
                    " effectively highlighting the exon structures."
        },
        {
            "key": "exon_scale",
            "annotation": "float",
            "default": "1",
            "note": "The degree of exon scaling. See intron_scale."
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


if __name__ == '__main__':
    pass

# trackplot Python chain API

Use the API instead of the CLI when you want reproducibility, no shell escaping, or to build
figures inside notebooks/scripts. API mirrors the CLI one-to-one.

## Table of contents
- [Quick start (chain)]#quick-start-chain
- [Class: Plot](#class-plot)
- [Method signatures](#method-signatures)
- [Single-cell example](#single-cell-example)

## Quick start (chain)
```python
from trackplot.plot import Plot

plot = Plot(logfile=None, backend="agg", font_family=None)
plot.set_annotation("example/example.sorted.gtf.gz", add_domain=False, show_gene=True) \
    .set_region("chr1", 1270656, 1284730, "+") \
    .add_density(path="example/bams/1.bam", category="bam", color="blue", show_site_plot=False) \
    .add_density(path="example/bws/2.bw", category="bw", color="green") \
    .add_line(path="example/bams/2.bam", category="bam", group="2", color="red", line_attrs={"linestyle": "dashed"}) \
    .add_heatmap(path="example/bams/1.bam", category="bam", group="1") \
    .add_igv(path="example/bams/3.bam", category="igv", label="igv", features={}) \
    .add_sites(1271656).add_sites(1271656).add_sites(1272656) \
    .add_focus("1272656-1272656:1275656-1277656") \
    .add_stroke("1275656-1277656:1277856-1278656@blue") \
    .plot("test.png", width=6, height=1, raster=True, dpi=300)
```
Each `add_*` method returns `self`, so calls chain.

## Class: Plot
```python
Plot(logfile=None, backend="agg", font_family=None)
```
- `logfile`: save progress logs to a file (default stdout).
- `backend`: matplotlib backend. `Agg`/`PDF` may drop protein domains; `Cairo` needs `cairocffi`.
- `font_family`: override font family (default system).

## Method signatures

### set_region(chromosome, start, end, strand="+")
Target locus. Returns `Plot`.

### set_annotation(gtf, add_domain=False, local_domain=False, domain_include=False,
domain_exclude=False, interval=None, interval_label=None, transcripts=None,
remove_empty_transcripts=False, choose_primary=False, color="black",
font_size=5, show_gene=False, show_id=False, exon_width=.3, show_exon_id=False, theme="blank")
- `theme`: one of `blank`, `ticks`, `ticks_blank`.
- `interval`/`interval_label`: add a custom BED feature track to the annotation.
- `transcripts`: list of transcript names/ids to draw; `choose_primary`: plot primary transcript only.

### add_density(path, category="bam", size_factor=None, label="", title="", barcode="",
barcode_groups=None, barcode_tag="BC", umi_tag="UB", library="fru",
density_by_strand=False, color="blue", font_size=8, show_junction_number=True,
junction_number_font_size=5, n_y_ticks=4, show_y_label=True, y_label="",
theme="ticks_blank", log_trans=None, show_site_plot=False, strand_choice=None,
only_customized_junction=False)
- `category`: `bam`/`bw`/`depth`/`bgz`.
- `library`: `fru`=fr-unstrand, `frf`=fr-firststrand, `frs`=fr-secondstrand.
- `log_trans`: `0`=none, `2`=log2, `10`=log10.
- `show_site_plot`: draw read-start density per strand (needs `--site-strand` equivalent via strand_choice).

### add_line(path, group="", category="bam", label="", title="", color=None,
line_attrs={}, **common)  — same file-loading kwargs as add_density.

### add_heatmap(path, group="", category="bam", size_factor=None, label="", title="",
color="viridis", font_size=8, show_y_label=True, theme="ticks_blank", do_scale=False,
clustering=False, clustering_method="ward", distance_metric="euclidean",
show_row_names=False, vmin=None, vmax=None, log_trans=None)
- `clustering_method` / `distance_metric` follow scipy (single/complete/average/weighted/centroid/median/ward; euclidean/cosine/correlation/...).

### add_igv(path, category="igv", label="", color=None, features={})
- `features`: read tags to overlay, e.g. `{"m6a": "ma", "polya": "pa", "real_strand": "rs"}`.

### add_hic(path, category="hic", label="", color="RdYlBu_r", transform=None, depth=None)

### add_stroke(text) / add_focus(text) / add_link(text)
- `stroke`: `start1-end1:start2-end2@color-label` (bottom stroke line; default red).
- `focus`: `100-200:300-400` (highlight region).
- `link`: `start1-end1:start1-end1@color` (bottom link between two sites).
- Also accept keyword form: `add_stroke(start=.., end=.., color="green", label="..")`.

### add_sites(*points) — comma/positional indicator lines; repeated points get a distinct color.

### plot(output, width=10, height=1, dpi=300, raster=False, **kwargs)
Render and save. `output`: png/pdf/svg/jpg. `raster=True` rasterizes heatmap/sites (smaller/faster PDF/SVG).

## Single-cell example
```python
plot = Plot()
plot.set_annotation("example/example.sorted.gtf.gz") \
    .set_region("chr1", 1270656, 1284730, "+") \
    .add_density(path="example/bams/sc.bam", category="bam",
                 barcode_groups={"AT2": {"AAACCTGCACCTCGTT-1"}},
                 barcode_tag="CB", umi_tag="UB", density_by_strand=True) \
    .plot("sc.pdf", height=1)
```
Equivalently pass a barcode list file via CLI `--barcode <file> --group-by-cell`.

---
name: trackplot
description: "Generate sashimi-style genome visualization plots (coverage, line, heatmap, IGV read-by-read, HiC, circRNA, motif) from BAM/bigWig/depth/HiC inputs. Use when the user wants to plot NGS data over a genomic region, make sashimi or intron-shrinkage plots, strand density, single-cell barcode-split density, protein domain tracks, or publish-ready PDF/PNG/SVG figures for a locus. Covers CLI usage, config TSV file formats, Docker, and the Python Plot chain API. Use when the user says trackplot, sashimi plot, coverage plot for a locus, or chr:start-end:strand plotting with bigWig or BAM."
license: BSD-3-Clause. LICENSE.txt has full terms.
---

# trackplot

trackplot is a pure-Python (>=3.8) sashimi-plot / locus-visualization framework. It draws
coverage, line, heatmap, individual-read (IGV), HiC, circRNA and motif tracks for a single
genomic region, and emits journal-ready PDF/PNG/SVG. Input is given as tab-separated config
files; output is one figure where each track maps to one config file.

Repo: https://github.com/ygidtu/trackplot · Docs: https://trackplot.readthedocs.io · DOI: 10.1371/journal.pcbi.1011477

## When to choose which track

| User wants | Use |
|-----------|-----|
| Coverage/sashimi of BAM or bigWig over a locus (junctions shown) | `--density` |
| Multi-sample coverage as lines (time course, conditions) | `--line` |
| Several samples side-by-side coverage blocks | `--heatmap` |
| Individual aligned reads, incl. long-read m6A/polyA marks | `--igv` |
| 2D contact matrix / HiC | `--hic` |
| Gene model, exon/intron, protein domains, custom beds | `-r` annotation + `--domain` / `--interval` |
| Circular RNA / back-splice highlight | `--density` + `--stroke` + `--link` |

For config-file formats (exact columns) and the Python chain API, see
[references/config_files.md](references/config_files.md) and
[references/python_api.md](references/python_api.md).

## Installation

Fastest (Linux/macOS x86 with glibc):
```bash
pip install trackplot                 # bigWig/bigBed/HiC support is optional, see below
```

Optional extras (enable formats that otherwise error out):
```bash
pip install pybigwig hicmatrix        # bigWig, bigBed, and .hic / .h5
```

Other supported installs: bioconda (`conda install -c bioconda -c conda-forge trackplot`),
source (`pip install -e .`), `uv` (`uv sync`), AppImage (Linux/WSL x86_64 only), Docker.

**Platform caveats** (from upstream docs):
- Windows, Apple-Silicon macOS and other ARM hosts often **cannot** install via PyPI/conda
  because of pysam/pybigwig/hicmatrix wheels. **Use the Docker image there.**
- If multi-processing causes `segment fault`, rerun with `-p 1` or use Docker.
- If you see `Please install pyBigWig and hicmatrix`, install the optional extras above.

```bash
# Docker (recommended on macOS/ARM/Windows)
docker pull ygidtu/trackplot
docker run --rm -v $PWD:$PWD -w $PWD ygidtu/trackplot --help
```

## Region format (required)

Every plot targets one region with `-e`:
```
chromosome_id:start:end:strand      e.g.  chr1:1270656-1284730:+
```
Strand is `+` or `-`. For strand-aware density use `--density-by-strand`.

## Common CLI workflow (density / sashimi plot)

```bash
trackplot \
  -e chr1:1270656-1284730:+ \
  -r example/example.sorted.gtf.gz \
  --density example/density_list.tsv \
  --show-junction-num \
  -o figure.pdf \
  --dpi 300 --width 10 --height 1 \
  -p 4
```

- `-r/--annotation`: GTF/GFF (both `transcript` and `exon` tags must be present).
  Sorted + bgzipped+tabix is fine but not required.
- `--density` / `--line` / `--heatmap` / `--igv` / `--hic`: each takes a config TSV path.
  Add as many different track types as needed in one invocation.
- `-o/--output`: `pdf`, `png`, `svg`, `jpg` all supported. Journal PDF: use vector; for
  large heatmaps/sites pass `--raster` to shrink file size / speed up rendering.

### Frequently used options
- `--show-junction-num` / `--show-mean-junction-num`: annotate intron junction read counts.
- `-t/--threshold`: drop low-abundance junctions (min count).
- `--show-site`: draw read-start position (site) marks on the density.
- `--focus 100-200:300-400`: highlight a region; `--stroke a-b:c-d@color-label`: bottom stroke line;
  `--link a-b:c-d@color`: bottom link between two sites; `--sites 12,34,56`: comma-separated indicator lines.
- `--intron-scale 0.5` / `--exon-scale 1`: shrink/expand introns (fixed introns: scale > 1).
- `--domain`: add protein domain track from UniProt/Ensembl (needs network) or `--local-domain <folder>` (UCSC bigBed).
- `--interval <bed.tsv>`: add custom feature track to the annotation.
- `--log [0|2|10|zscore]`: log-transform the y axis. `--normalize-format [count|cpm|rpkm]`: normalize BAM.
- `--color-factor N`: color by a categorical column of the config file (`LUAD|red` → label LUAD, color red).
- `--width/--height/--dpi/--backend/--font-size/--title/--font`: output/figure styling.

### Single-cell BAM (barcode-split density/line)
Requires a barcode list and tags (10x default: `--barcode-tag CB --umi-tag UB`):
```bash
trackplot \
  -e chr1:1270656-1284730:+ \
  -r example/example.sorted.gtf.gz \
  --density example/density_list.tsv \
  --barcode example/barcode_list.tsv \
  --group-by-cell \
  -o sc.pdf
```
Barcode list columns: `bam  barcode  cell_type(optional)  color(optional)`.

## Config-file formats (tab-separated)

Header line starts with `#`; comment/blank lines ignored. Columns:

```
# density     # filepath  category   label(optional)  color(optional)
# line        # filepath  category   group(optional)  color(optional)
# heatmap     # filepath  category   group(optional)  color(optional)
# igv         # filepath  category   label(optional)  color(optional)
# hic         # filepath  category   label(optional)  color(optional)  transform(optional)  depth(optional)  domain(optional)
# interval    # file_location  label
# custom-junction  junctions  <bam-or-aliases...>   then   <junction-id>  <count-per-column...>
```
- `category` is one of `bam`, `bw`, `bed`, `depth`, `hic`, `igv`, `bed3/6/12`, etc.
- For `bam` in density/line/heatmap you may append library / total-read columns; see
  [references/config_files.md](references/config_files.md) for the full column table.

## Docker usage notes
- **Absolute paths required** inside the container. Convert relative config paths:
  ```bash
  grep -v '^#' example/density_list.tsv | while read l; do echo "$PWD/${l}"; done > abspath.tsv
  docker run -v $PWD:$PWD -w $PWD --rm ygidtu/trackplot -e chr1:...:... --density abspath.tsv -o out.pdf
  ```
- Mount your data dir with `-v` and set `-w` to it so paths match.

## Web UI (optional)
Start a local server for a browser-based plot builder:
```bash
trackplot --start-server --host 127.0.0.1 --port 5000 --plots ./plots     # --plots required for AppImage
# docker: docker run -v $PWD/example:/data -v $PWD/plots:/plots -p 5000:5000 ygidtu/trackplot --start-server --data /data --plots /plots
```
Region must match `chromosome_id:start_site-end_site:strand`.

## Troubleshooting
- `#REF!`-like import errors / pysam build issues → use Docker (see platform caveats).
- Missing bigWig/HiC support → `pip install pybigwig hicmatrix`.
- Protein domains missing in output → avoid `Cairo` backend for `--domain`; use `Agg`/`PDF`.
- Slow PDF/SVG with heatmaps/sites → add `--raster`.
- Want reproducibility / no shell escaping → prefer the Python chain API; see
  [references/python_api.md](references/python_api.md).

## Do not
- Do not pass a region without strand; do not omit `-r` when junctions/annotation is expected.
- Do not hardcode coverage values in Python then plot them — let trackplot read BAM/bigWig directly.
- Do not expect `vite build`-style type checks here; this is a plotting tool, not a web framework.

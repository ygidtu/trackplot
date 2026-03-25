# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Trackplot is a Python tool for visualizing next-generation sequencing (NGS) data (DNA-seq, RNA-seq, single-cell RNA-seq, full-length sequencing). It generates publication-quality figures using matplotlib.

## Development Commands

### Setup
```bash
# Install in development mode
pip install -e .

# Install with all dependencies
pip install -r requirements.txt
pip install -e .
```

### Running the Tool
```bash
# CLI mode
python main.py --help

# Or after pip install
trackplot --help

# Example basic plot
trackplot -e chr1:1270656-1284730:+ -r example/example.sorted.gtf.gz -o output.png
```

### Web Interface
```bash
# Start web server
trackplot --start-server --host 0.0.0.0 --port 5000 --plots ./plots

# Build frontend (requires pnpm)
cd web
pnpm install
pnpm build
```

### Docker
```bash
# Build Docker image
docker build -t trackplot .

# Run with Docker
docker run --rm -v $PWD:/data trackplot -e chr1:1270656-1284730:+ -r /data/example/example.sorted.gtf.gz -o /data/output.png
```

## Architecture Overview

### Core Module Structure (`trackplot/`)

**Entry Points:**
- `main.py` - CLI entry point, calls `cli.py`
- `cli.py` - Command-line argument parsing and file list processing
- `server.py` - Flask web server for the UI

**Core Plotting:**
- `plot.py` - Main `Plot` class that orchestrates figure creation; contains `PlotInfo` for tracking plot datasets
- `plot_func.py` - Individual plotting functions (`plot_density`, `plot_reference`, `plot_heatmap`, etc.)

**Data Model (`base/`):**
- `GenomicLoci.py` - Core genomic coordinate class (chrom, start, end, strand); most classes inherit from this
- `Transcript.py` - Transcript annotation with exons
- `Protein.py` - Protein domain handling
- `ReadDepth.py` - Coverage depth data container
- `Junction.py` - Splice junction representation
- `CoordinateMap.py` - Coordinate transformation utilities
- `Readder.py` - File reading utilities
- `Stroke.py` - SVG-like stroke annotations

**File Format Handlers (`file/`):**
- `File.py` - Base class for all input files with `load()` interface
- `Annotation.py` - GTF/GFF annotation parser
- `Bam.py` - BAM alignment file handling
- `Bigwig.py` - BigWig coverage files
- `BedGraph.py` - BedGraph format
- `Depth.py` - Samtools depth output
- `ReadSegments.py` - Individual read visualization
- `HiCMatrixTrack.py` - Hi-C interaction matrices
- `ATAC.py` - ATAC-seq data
- `Fasta.py` - Fasta sequence files
- `Motif.py` - Motif/sequence annotations
- `Junction.py` - Custom junction file parsing

**Configuration (`conf/`):**
- `config.py` - Constants (colors, normalization methods, etc.)
- `ui.py` - Web UI parameter definitions
- `drawing.py` - Drawing configuration classes
- `DomainSetting.py` - Protein domain display settings

**Annotation (`anno/`):**
- `theme.py` - Matplotlib theme/styling
- `AxLabel.py` - Axis labeling utilities

### Web Interface (`web/`)

Vue 3 + TypeScript + Vite frontend:
- Built with Element Plus UI components
- Communicates with Flask backend in `server.py`
- Build output goes to `ui/` (generated, not in repo)

### Key Design Patterns

1. **File Loading Pattern**: All file handlers inherit from `File` base class and implement `load(region, **kwargs)` method
2. **GenomicLoci as Foundation**: Most data structures extend `GenomicLoci` for chromosome/start/end/strand handling
3. **Plot Composition**: `Plot` class aggregates `PlotInfo` objects, each representing one track
4. **Multiprocessing**: Data loading uses `multiprocessing.Pool` for parallelization (see `__load__` function in `plot.py`)

## File Formats

The tool supports multiple input formats:
- **Density tracks**: BAM, BigWig, BedGraph, Depth files
- **Annotation**: GTF, GFF, BED
- **Read visualization**: BAM (read-by-read mode)
- **Heatmaps**: Hi-C matrices (hicmatrix format)
- **Intervals**: BED files for marking regions

## Important Notes

- Python >=3.11 required (defined in `pyproject.toml`)
- Uses `loguru` for logging
- Heavy dependencies: `pysam`, `pybigwig`, `hicmatrix`, `matplotlib`, `scipy`
- Windows/ARM users should use Docker due to pysam/pybigwig compilation issues
- If segment faults occur during multiprocessing, use `-p 1` flag

## CLI Optimizations

### Lazy Import Optimization (64% improvement)

The original `cli.py` imported `server` module at the top level, causing Flask to load even in pure CLI mode. Now merged into `cli.py`:
```python
# Import only when needed (inside main())
if kwargs["start_server"]:
    from trackplot.server import run
    run(...)
```

**Benchmark Results:**
- Import time: 3.432s → 1.238s (64% faster)
- Server module loaded: True → False (not loaded in CLI mode)
- Flask loaded: True → False (not loaded in CLI mode)

### process_file_list Refactoring

The 200+ line nested if-elif was refactored with centralized validation (`_CATEGORY_VALIDATORS` dict) and cleaner variable naming (`cat` instead of shadowing `category`). Use `category == "x"` instead of `category in ["x"]` for single-value checks.

### Remaining Optimization Opportunities

1. **`os.path.exists()` bottleneck**: `FileList.__init__` calls `os.path.exists()` on every entry. For large file lists (1000+), this dominates runtime (~3ms for 1000 entries). Could defer validation to load time.

2. **Repeated `groups.get()` pattern in line parsing**: Each branch repeats `groups[line[2]] = groups.get(line[2], 0) + 1`. Could extract once before the if-elif chain.

3. **Import graph**: `cli.py` still imports `ATAC` at top level even when not used. Could lazy-import similar to server.

4. **`main()` function is ~300 lines**: The click-decorated `main()` function handles all plotting dispatch. Extracting per-category handlers (e.g., `_handle_density()`, `_handle_heatmap()`) would improve readability.

### Benchmark Scripts

```bash
# Test import time improvement
python bench_opt.py

# Test parser performance (compares cli.py vs cli_optimized.py)
python bench_parser.py

# Test parsing logic only
python bench_micro.py
```

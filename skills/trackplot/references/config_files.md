# trackplot config-file formats

All input track files are tab-separated. Lines beginning with `#` and blank lines are ignored.
The first column is always the **file path** (relative to CWD, or absolute when using Docker).
`category` names: `bam`, `bw` (bigWig), `bgz` (samtools depth, bgzipped+tabix), `bed`, `hic`/`h5`,
and IGV read formats `bed3/bed6/bed12`.

## Table of contents
- [density (sashimi / coverage)]#density
- [line](#line)
- [heatmap](#heatmap)
- [igv (read-by-read)]#igv
- [hic](#hic)
- [interval (annotation features)]#interval
- [customized-junction](#customized-junction)
- [junction counts](#junction-counts)
- [barcode (single-cell)]#barcode
- [motif](#motif)
- [colorFactor / color conventions](#colorfactor--color-conventions)

## density
```
# filepath  category   label(optional)  color(optional)
example/bams/1.bam  bam
example/bams/2.bam  bam  2bam
example/bams/3.bam  bam  3bam  blue
example/bws/2.bw    bw   bw green
```
- Accepts `bam`, `bw`, or `bgz` depth files.
- Extra columns for **single-cell bam** (5th/6th): `library` / `total reads` — see [barcode](#barcode).

## line
```
# filepath  category   group(optional)  color(optional)
example/bams/1.bam  bam  bam
example/bams/2.bam  bam  bam
```
Same category semantics as density; draws each file as a line, grouped by the 3rd column.

## heatmap
```
# filepath  category   group(optional)  color(optional)
example/bams/1.bam  bam  bam
example/bws/1.bw    bw   bw  YlOrBr
```
Files in the same `group` render in the same heatmap panel.

## igv
```
# filepath  category   label(optional)  color(optional)
example/SRX9697989.corrected_reads.bed.gz  igv  bed12  blue
example/bams/0.bam                           igv  bam
```
- `category` selects read format (`bed12`, `bed6`, `bed3`, `bam`).
- Long-read marks via BAM tags on the command line:
  `--m6a ma` (tag `ma:i`, modification site), `--polya pa` (tag `pa:f`, poly(A) length),
  `--rs rs` (real strand tag, required to draw poly(A)).

## hic
```
# filepath  category   label  color       transform  depth   domain
example/Li_et_al_2015.h5  hic  Li_hic  RdYlBu_r  2          50000   example/..._domains.bed.gz
```
- `transform`: `0`=none, `2`=log2, `10`=log10, `e`=ln (same as `--log`).
- `depth`: larger value → taller y-axis for that track.
- `domain` (optional): TAD bed file from `hicFindTADs`.

## interval
Added to the annotation track with `--interval`:
```
# file_location                label
example/PolyASite.chr1...simple.bed.gz  polyAS
```

## customized-junction
Draws user-supplied junctions (counts per bam/alias):
```
junctions  2bam  3bam
chr1:1271656-1272656  100  200
chr1:1273656-1274656  100  200
```
First row is the header (junction id + one column per input). Each following row: a junction id
(`chr:start-end`) and its read count in each sample.

## junction counts
Used with `--show-junction-num`; first column is the input file name, then counts:
```
junctions  1.bam
chr1:1270655-1271655:+  500
chr1:1277757-1277857:+  501
```

## barcode
For single-cell BAM density/line. Columns:
```
# bam  barcode  cell_type(optional)  cell_color(optional)
sc  AAACCTGCACCTCGTT-1  AT2  #A6DCC2
```
Pair with `--barcode <file> --barcode-tag CB --umi-tag UB --group-by-cell`.
Default tags are 10x Genomics (`CB` = cell barcode, `UB` = UMI).

## motif
Customized bedGraph; first three columns chrom/start/end, next four are ATCG weights:
```
# chromosome  start  end  A_weight  T_weight  C_weight  G_weight
chr1  100  101  0.1  0.2  -0.3  -0.4
```
bgzip + tabix it, then `--motif <file> --motif-region 1270756-1270760`.

## colorFactor / color conventions
- In any list file, put a label in the appropriate column and set `--color-factor <col>` (1-based)
  so trackplot colors by that category: `path LUAD` + `--color-factor 2` → LUAD gets an auto color.
- Inline override with `name|hexcolor`: `path LUAD|red` → label LUAD, color red; `path LUSC|#000000`.

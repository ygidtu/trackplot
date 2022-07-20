
# Interactive Usage

## Example Usage

```python
from sashimi.plot import Plot
```

The interactive API could be used as normal functions or calling as process chain.

1. init plot


```python
plot = Plot()
```

2. set reference


```python
plot.set_reference(
    "../example/example.sorted.gtf.gz",
    add_domain=True,
    interval="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.bed.gz",
    interval_label="polyA",
    show_gene=True,
    color="pink",
)
```

    building of index for ../example/example.sorted.gtf.gz failed
    Guess gtf needs to be sorted





    <sashimi.plot.Plot at 0x7f91e1074100>



3. setup plotting parameters


```python
plot.set_reference(
    "../example/example.sorted.gtf.gz",
    add_domain=True,
    interval="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.bed.gz",
    interval_label="polyA",
    show_gene=True,
    color="pink",
).add_interval(
    interval="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.simple.bed.gz",
    interval_label="polyAS"
).set_region(
    "chr1", 1270656, 1284730, "+"
).add_density(
    path="../example/bams/1.bam",
    category="bam",
    color="blue",
    show_side_plot=True,
).add_density(
    path="../example/bws/2.bw",
    category="bw",
    color="green"
).add_line(
    path="../example/bams/1.bam",
    category="bam",
    group="1",
    line_attrs={"lw": 3}
).add_line(
    path="../example/bams/2.bam",
    category="bam",
    group="2",
    color="red",
    line_attrs={"linestyle": "dashed"}
).add_heatmap(
    path="../example/bams/1.bam",
    category="bam",
    group="1",
).add_heatmap(
    path="../example/bams/2.bam",
    category="bam",
    group="1"
).add_igv(
    path="../example/bams/3.bam",
    features={
        "m6a": "ma",
        "real_strand": "rs",
        "polya": "pa"
    },
    category="igv",
    label="igv"
).add_igv(
    path="../example/SRX9697989.corrected_reads.bed.gz",
    category="igv",
    label="bed12"
).add_sites(
    1270656 + 1000
).add_sites(
    1270656 + 1000
).add_sites(
    1270656 + 2000
).add_focus(
    f"{1270656 + 2000}-{1270656 + 3000}"
).add_focus(
    f"{1270656 + 5000}-{1270656 + 7000}"
).add_stroke(
    f"{1270656 + 5000}-{1270656 + 7000}:{1270656 + 7200}-{1270656 + 8000}@blue"
).add_stroke(
    start=1270656 + 7500,
    end=1270656 + 8200,
    color="green",
    label="test"
)
```

4. save figures


```python
plot.plot("test_plot.png", fig_width=6, fig_height=2, raster=True)
```

## API documentation

```python
def set_region(self, chromosome: str, start: int, end: int, strand: str = "+")
```

change the plot region

- chromosome:
- start:
- end:
- strand:

Return

---

```python
def set_reference(self, gtf: str,
                  add_domain: bool = False,
                  local_domain: Optional[str] = False,
                  interval: Optional[str] = None,
                  interval_label: Optional[str] = None,
                  transcripts: Optional[List[str]] = None,
                  remove_empty_transcripts: bool = False,
                  color: Optional[str] = None,

                  # transcripts related parameters
                  font_size: int = 5,
                  show_gene: bool = False,
                  show_id: bool = False,
                  reverse_minus: bool = False,
                  exon_width: float = .3,
                  show_exon_id: bool = False,
                  theme: str = "blank"
                  )
```

add transcripts to this region

- gtf: path to gtf file
- add_domain:
- local_domain:
- interval:
- interval_label:
- font_size: the size of transcript id, name
- transcripts: the list of name or ids of transcripts to draw
- remove_empty_transcripts: whether to remove transcripts without any exons
- color: the color of exons
- show_gene: whether to show gene name/id
- show_id: show gene id or gene name
- reverse_minus: whether to remove strand of transcripts
- theme: the theme of transcript
- exon_width: the height of exons
- show_exon_id: whether to show exon id

Returns `Plot`

---

```python
 def add_density(self,
                path: str,
                category: str,

                # file loading parameters
                label: Union[str, List[str]] = "",
                title: str = "",
                barcodes: Optional[Set[str]] = None,
                barcode_tag: str = "BC",
                umi_tag: str = "UB",
                library: str = "fr-unstrand",

                # plotting parameters
                color="blue",
                font_size: int = 8,
                show_junction_number: bool = True,
                junction_number_font_size: int = 5,
                n_y_ticks: int = 4,
                distance_between_label_axis: float = .1,
                show_y_label: bool = True,
                y_label: str = "",
                theme: str = "ticks_blank",

                # side plot parameters
                show_side_plot: bool = False,
                strand_choice: Optional[str] = None,
                )
```

add density object to plot

- path: the path to input file
- category: the input file type
- show_side_plot: draw the density distribution of reads from different strand
- label: the label of input file
- title: the title of input file
- barcodes: list of required barcodes
- barcode_tag: cell barcode tag
- umi_tag: umi barcode tag
- library: fr-unstrand
- font_size: the font size for ticks, y-axis label and title
- show_junction_number: whether to show the number of junctions
- distance_between_label_axis: distance between y-axis label and y-axis ticks
- n_y_ticks: number of y ticks
- junction_number_font_size:
- color: color for this density plot
- show_y_label: whether to show y-axis label
- y_label: the text of y-axis title
- theme: the theme name
- strand_choice: the strand to draw on side plot

Returns `Plot`

---

```python
def add_heatmap(self,
                path: str,
                group: str,
                category: str = "bam",

                # file loading parameters
                label: Union[str, List[str]] = "",
                title: str = "",
                barcodes: Optional[Set[str]] = None,
                barcode_tag: str = "BC",
                umi_tag: str = "UB",
                library: str = "fr-unstrand",

                # plotting parameters
                color="viridis",
                font_size: int = 8,
                distance_between_label_axis: float = .1,
                show_y_label: bool = True,
                theme: str = "ticks",
                do_scale: bool = False,
                clustering: bool = False,
                clustering_method: str = "ward",
                distance_metric: str = "euclidean",
                )
```

add multiple objects for a group of heatmap

- path: path to input files
- group: the heatmap group
- category: file category corresponding to input file
- label: the label of input file
- title: the title of input file
- barcodes: list of required barcodes
- barcode_tag: cell barcode tag
- umi_tag: umi barcode tag
- library: fr-unstrand
- color: color for this density plot
- show_y_label: whether to show y-axis label
- theme: the theme name
- font_size:
- distance_between_label_axis:
- do_scale: whether to scale the matrix
- clustering: whether reorder matrix by clustering
- clustering_method: same as  scipy.cluster.hierarchy.linkage
- distance_metric: same as scipy.spatial.distance.pdist
- color: used for seaborn.heatmap, see: https://matplotlib.org/3.5.1/tutorials/colors/colormaps.html
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
            'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
            'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'

Returns `Plot`


---

```python
def add_line(self,
             path: str,
             group: str,
             category: str = "bam",

             # file loading parameters
             label: Union[str, List[str]] = "",
             title: str = "",
             barcodes: Optional[Set[str]] = None,
             barcode_tag: str = "BC",
             umi_tag: str = "UB",
             library: str = "fr-unstrand",

             # plotting parameters
             color="blue",
             font_size: int = 8,
             distance_between_label_axis: float = .1,
             show_y_label: bool = True,
             line_attrs: Optional[Dict] = None,
             theme: str = "ticks_blank",
             n_y_ticks: int = 4,
             show_legend: bool = False,
             legend_position: str = "upper right",
             legend_ncol: int = 0
             )
```

add multiple objects for a group of heatmap

- path: path to input files
- group: the heatmap group
- category: file category corresponding to input file
- label: the label of input file
- title: the title of input file
- barcodes: list of required barcodes
- barcode_tag: cell barcode tag
- umi_tag: umi barcode tag
- library: fr-unstrand
- distance_between_label_axis: distance between y-axis label and y-axis ticks
- n_y_ticks: number of y ticks
- color: color for this density plot
- show_y_label: whether to show y-axis label
- theme: the theme name
- font_size: font size in this plot
- line_attrs: the additional attributes to control the line, usd by matpltolib.axes.Axes.plot
- show_legend: whether to show legend
- legend_position: the position of legend
- legend_ncol: the number of columns in legend

Returns `Plot`

---

```python
def add_igv(
            self,
            path: str,
            category: str = "igv",
            label: str = "",
            exon_focus: Optional[str] = None,

            # file loading parameters
            library: str = "fr-unstrand",
            features: Optional[dict] = None,
            deletion_ignore: Optional[int] = True,
            del_ratio_ignore: float = .5,

            # plotting parameters
            exon_color: Optional[str] = None,
            intron_color: Optional[str] = None,
            feature_color: Optional[str] = None,
            exon_width: float = .3,
            font_size: int = 8,
            n_y_ticks: int = 1,
            distance_between_label_axis: float = .1,
            show_y_label: bool = True,
            theme: str = "ticks_blank"
)
```

Add igv-like plot into track

- path: path to input files
- category: file category for the input file
- library: fr-unstrand
- features:
- exon_focus:
- deletion_ignore:
- del_ratio_ignore:
- label:
- exon_color:
- intron_color:
- feature_color:
- exon_width:
- font_size:
- n_y_ticks:
- distance_between_label_axis:
- show_y_label:
- theme:

Returns `Plot`


---

```python
def add_sites(self, sites)
```

highlight specific sites

- sites: string in 100,200 format or int

Returns `Plot`

---

```python
def add_focus(self, focus: Optional[str], start: int = 0, end: int = 0)
```

set focus region
- focus: string in 100-200:300-400
- start: start site
- end: end site

Returns `Plot`

--- 

```python
def add_stroke(
                self,
                stroke: Optional[str] = None,
                start: int = 0,
                end: int = 0,
                label: str = "",
                color: str = "black"
)
```

- stroke: string format of stroke,  eg: 100-200:200-300@blue
- start: start position
- end: position
- label: stroke label
- color: the color of stroke

Returns `Plot`

---

```python
def set_sequence(self, fasta: str)
```

set sequence info
- fasta: path to indexed fasta file

Returns `Plot`


---

```python
def add_interval(self, interval: str, interval_label: str)
```

---

```python
def plot(self,
         output: Optional[str] = None,
         reference_scale: Union[int, float] = .25,
         stroke_scale: Union[int, float] = .25,
         dpi: int = 300,
         fig_width: Union[int, float] = 0,
         fig_height: Union[int, float] = 0,
         raster: bool = False,
         *args, **kwargs)
```
save image

- output: if output is empty then show this image by plt.showfig
- reference_scale: to adjust the size of reference plot
- stroke_scale: to adjust the size of stroke plot
- dpi: the dpi of saved plot
- fig_width: the width of figure, if width == 0, the let matplotlib decide the size of image
- fig_height: the height of figure, if height == 0, the let matplotlib decide the size of image
- raster: plot rasterizer side plot


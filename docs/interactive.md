
# Interactive Usage

## Example Usage

The interactive API could be used as normal functions or calling as process chain.

###### 1. init plot

```python
from sashimi.plot import Plot

plot = Plot()
```

###### 2. set reference

```python
plot.set_reference(
    "../example/example.sorted.gtf.gz",                                             # path to gtf file
    add_domain=True,                                                                # whether add domain information into reference track
    interval="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.bed.gz",       # path to list of interval files in bed format, 1st column is path to file, 2nd column is the label
    interval_label="polyA",                                                         # the label of added interval
    show_gene=True,                                                                 # show gene id
    color="pink",                                                                   # the color of exons
)
```

###### 3. setup plotting parameters

```python
plot.add_interval(
    interval="../example/PolyASite.chr1.atlas.clusters.2.0.GRCh38.96.simple.bed.gz",          # path to list of interval files in bed format, 1st column is path to file, 2nd column is the label
    interval_label="polyAS"                                                                   # the label of added interval
).set_region(
    chromosome="chr1", start=1270656, end=1284730, strand="+"                                 # chromosome, start site, end site and strand were required
).add_density(
    path="../example/bams/1.bam",                                                             # path to input file
    category="bam",                                                                           # the category of given file
    color="blue",                                                                             # color of this density
    show_site_plot=True,                                                                      # whether to show site plot
).add_density(                                                                                # another density plot
    path="../example/bws/2.bw",
    category="bw",
    color="green"
).add_line(
    path="../example/bams/1.bam",                                                             # path to input file
    category="bam",                                                                           # the category of given file
    group="1",                                                                                # the group of this line, used to control color etc.
    color="blue",                                                                             # color of this line and other lines belong to group 1
    line_attrs={"lw": 3}                                                                      # additional parameters in dict to control the layout of line
).add_line(                                                                                   # another line
    path="../example/bams/2.bam",
    category="bam",
    group="2",
    color="red",
    line_attrs={"linestyle": "dashed"}                                                        # additional parameters in dict to control the layout of line, for instance this changes the line style
).add_heatmap(
    path="../example/bams/1.bam",                                                             # path to input file
    category="bam",                                                                           # the category of given file
    group="1",                                                                                # the group of this file, only files belong to same group will be drawn in same heatmap.
).add_heatmap(                                                                                # another file to heatmap 1
    path="../example/bams/2.bam",
    category="bam",
    group="1"
).add_igv(
    path="../example/bams/3.bam",                                                             # path to input file
    features={
        "m6a": "ma",
        "real_strand": "rs",
        "polya": "pa"
    },
    category="igv",                                                                           # the category of given file
    label="igv"                                                                               # the label of this plot
).add_igv(
    path="../example/SRX9697989.corrected_reads.bed.gz",
    category="igv",
    label="bed12"
).add_sites(
    1270656 + 1000                                                                            # the highlight site
).add_sites(
    1270656 + 1000                                                                            # the repeat highlight site will show in different color with normal highlight site
).add_sites(
    1270656 + 2000
).add_focus(
    f"{1270656 + 2000}-{1270656 + 3000}"                                                      # the focus region in start_site-end_site format
).add_focus(
    f"{1270656 + 5000}-{1270656 + 7000}"
).add_stroke(
    f"{1270656 + 5000}-{1270656 + 7000}:{1270656 + 7200}-{1270656 + 8000}@blue"               # the stroke in start_site-end_site:start_site-end_site@color format, this will add 2 strokes and the last on will be blue
).add_stroke(                                                                                 # add stroke with named parameters
    start=1270656 + 7500,                                                                     # the start site of added stroke
    end=1270656 + 8200,                                                                       # the end site of added stroke
    color="green",                                                                            # the color of added stroke
    label="test"                                                                              # the label of added stroke
)
```

###### 4. save figures

```python
plot.plot("test_plot.png", fig_width=6, fig_height=2, raster=True)
# save plot into test_plot.png, with given fig width and height, the raster=True, will reduce layers in pdf or svg format
```

## API documentation


### set_region

change the plot region

```python
def set_region(self, chromosome: str, start: int, end: int, strand: str = "+")
```

- chromosome: the chromosome of given region
- start: the start site of given region
- end: the end site of given region
- strand: the strand of given region [optional]

Return

---

### set_reference

add transcripts into track

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

- gtf: path to gtf file
- add_domain: whether to add domain
- local_domain: whether add domain information into reference track
- interval: path to list of interval files in bed format, 1st column is path to file, 2nd column is the label
- interval_label: the label of added interval
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

### add_density

add a density plot into track

```python
 def add_density(self,
                path: str,
                category: str,

                # file loading parameters
                label: Union[str, List[str]] = "",
                title: str = "",
                barcode: str = "",
                barcode_groups: Dict[str, Set[str]] = None,
                barcode_tag: str = "BC",
                umi_tag: str = "UB",
                library: str = "fr-unstrand",
                density_by_strand: bool = False,

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

                # site plot parameters
                show_site_plot: bool = False,
                strand_choice: Optional[str] = None,
                 
                only_customized_junction: bool = False
                )
```

- path: the path to input file
- category: the input file type
- show_site_plot: draw the density distribution of reads from different strand
- label: the label of input file
- title: the title of input file
- barcode: key of barcode barcode_groups
- barcode_groups: dict contains barcodes by groups; key -> Set[str]
- barcode_tag: cell barcode tag
- umi_tag: umi barcode tag
- library: fr-unstrand/fr-strand or fru/frs for short 
- density_by_strand: whether to draw density plot in strand-specific manner.
- font_size: the font size for ticks, y-axis label and title
- show_junction_number: whether to show the number of junctions
- distance_between_label_axis: distance between y-axis label and y-axis ticks
- n_y_ticks: number of y ticks
- junction_number_font_size:
- color: color for this density plot
- show_y_label: whether to show y-axis label
- y_label: the text of y-axis title
- theme: the theme name
- strand_choice: the strand to draw on site plot
- only_customized_junction: only work with bam files, only draw customized junctions

Returns `Plot`

---

### add_heatmap

add a heatmap based on a group of objects into track

```python
def add_heatmap(self,
                path: str,
                group: str,
                category: str = "bam",

                # file loading parameters
                label: Union[str, List[str]] = "",
                title: str = "",
                barcode: str = "",
                barcode_groups: Dict[str, Set[str]] = None,
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

- path: path to input files
- group: the heatmap group
- category: file category corresponding to input file
- label: the label of input file
- title: the title of input file
- barcode: key of barcode barcode_groups
- barcode_groups: dict contains barcodes by groups; key -> Set[str]
- barcode_tag: cell barcode tag
- umi_tag: umi barcode tag
- library: fr-unstrand/fr-strand or fru/frs for short 
- color: color for this density plot
- show_y_label: whether to show y-axis label
- theme: the theme name
- font_size: the font size in plot
- distance_between_label_axis: the distance between y-axis label and y-axis
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

### add_line

add a line plot based on a group of objects into track

```python
def add_line(self,
             path: str,
             group: str,
             category: str = "bam",

             # file loading parameters
             label: Union[str, List[str]] = "",
             title: str = "",
             barcode: str = "",
             barcode_groups: Dict[str, Set[str]] = None,
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

- path: path to input files
- group: the heatmap group
- category: file category corresponding to input file
- label: the label of input file
- title: the title of input file
- barcode: key of barcode barcode_groups
- barcode_groups: dict contains barcodes by groups; key -> Set[str]
- barcode_tag: cell barcode tag
- umi_tag: umi barcode tag
- library: fr-unstrand/fr-strand or fru/frs for short
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

### add_igv

add an igv-like plot into track

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

### add_sites

add multiple highlight sites into track

```python
def add_sites(self, sites)
```

highlight specific sites

- sites: string in 100,200 format or int

Returns `Plot`

---

### add_focus

add multiple highlight background into track

```python
def add_focus(self, focus: Optional[str], start: int = 0, end: int = 0)
```

set focus region
- focus: string in 100-200:300-400
- start: start site
- end: end site

Returns `Plot`

--- 

### add_stroke

add multiple highlight region under transcripts into track

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

### set_sequence

display the corresponding sequence under x-axis 

```python
def set_sequence(self, fasta: str)
```

set sequence info
- fasta: path to indexed fasta file

Returns `Plot`


---

### add_interval

```python
def add_interval(self, interval: str, interval_label: str)
```

- interval: path to local interval file
- interval_label: the label of given interval

Returns `Plot`

---

### add_manual

draw line or density plot based on manually added data.

```python
def add_manual(self, data: np.array, image_type: str = "line", 
               label: str = "", group: str = "", color: str = "blue", 
               font_size: int = 8, n_y_ticks: int = 1, 
               show_y_label: bool = True, theme: str = "ticks_blank",)
```

- data: the manual data object, should be 1D np.array with same length of target region
- image_type: the plotting type, one of `line` [default], `density`, `heatmap`
- label: the label of given manual data
- group: the group of given data, used for `line` or `heatmap` plot
- color: the color of given data
- font_size: the font size of this plot
- n_y_ticks: the number of y ticks to show
- show_y_label: whether to show y label
- theme: the theme of this plot


Returns `Plot`

---

### add_motif

draw motif based on bedGraph file

```python
def add_motif(self, path: str, category: str = "motif", motif_region: GenomicLoci = None,
              width: float = 0.8, theme: str = "blank")
```

- path: the path to tabix indexed bedGraph file, first 3 columns is chromosome, start and end site, the rest 4 columns is scores for ATCG.
- motif_region: to specify the position of motif
- width: the width of characters

Returns `Plot`

### merge_by_cell

```python
def merge_by_cell(self)
```

This is used to merge input files by label, for instance, 

```python
p = Plot()
p.add_density('bam', label="cell1")
p.add_density('bam', label="cell2")
p.add_density('bam1', label="cell1")

# After that, 3 plot obj saved in p.plots

p.merge_by_cell()
# the bam and bam1 with label == 'cell1' merged into one
```

---

### plot

save/show the final image

```python
def plot(self,
            output: Optional[str] = None,
            dpi: int = 300,
            width: Union[int, float] = 0,
            height: Union[int, float] = 0,
            raster: bool = False,
            same_y: bool=False,
            remove_duplicate_umi: bool = False,
            threshold: int = 10,
            included_junctions = ["chr1:1-100"],
            return_image: Optional[str] = None,
            
            n_jobs: int = 1,
            normalize_format: str = "count",

            intron_scale: float=.5,
            exon_scale: float = 1,
         
            reference_scale: float = .25,
            stroke_scale: float = .25,
            sc_height_ratio: Optional[Dict[str, float]] = None,
            distance_between_label_axis: float = .3,
        )
```

- output: if output is empty then show this image by plt.showfig
- dpi: the dpi of saved plot
- width: the width of figure, if width == 0, the let matplotlib decide the size of image
- height: the height of figure, if height == 0, the let matplotlib decide the size of image
- raster: plot rasterizer site plot
- same_y: whether the density plots share same y-axis range
- remove_duplicate_umi: drop duplicated UMIs by barcode
- threshold: threshold to filter low abundance junctions
- included_junctions: the list of junctions to draw, the junction should be `chrom:start-end` format string
- return_image: used for interactive ui, this parameter takes `png` or `pdf`, then will return corresponding format of image in bytes array
- n_jobs: the number of processes to use while loading data, recommended for huge number or size of input files
- normalize_format: used to normalized input data, should be one of `count`[default], `cpm` or `rpkm`, only worked for bam file
- intron_cale: used to control the plotting scale of introns, the introns only half size by default
- exon_cale: used to control the plotting scale of exons, the introns only half size by default
- reference_scale: to adjust the max size of reference plot, the references only occupy at most 1/4 of figure height by default
- stroke_scale: to adjust the max size of stroke plot, the references only occupy at most 1/4 of figure height by default
- sc_height_ratio: to adjust the relative height of single cell related plots, including single cell density and heatmap
- distance_between_label_axis: to adjust the distance between y-axis label and y-axis ticks
- n_jobs: load data in how many processes

Returns `None` or `io.BytesIO`
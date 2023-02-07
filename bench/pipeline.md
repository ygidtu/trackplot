# Manual Pipeline

## MISO
### Run miso pipeline before plotting

The miso pipeline must run in advance, due to build-in sashimi_plot required the miso output files.

```bash
conda activate misopy

# generate miso index
cd ref
zcat Homo_sapiens.GRCh38.101.chr.gff3 > Homo_sapiens.GRCh38.101.chr.gff3
index_gff --index Homo_sapiens.GRCh38.101.chr.gff3 miso_index/

mkdir miso
for i in $(ls STAR/*.Aligned.sortedByCoord.out.bam);
do
  echo $i
  mkdir miso/${$(basename $i)/.Aligned.sortedByCoord.out.bam/""}
  miso --run ref/miso_index/ $i --output-dir miso/${$(basename $i)/.Aligned.sortedByCoord.out.bam/""} --read-len 44 -p 20
done
```

The plot settings for `sashimi_plot` located in `miso/settings.txt`
```yaml
[data]
# directory where BAM files are
bam_prefix = ./STAR/
# directory where MISO output is
miso_prefix = ./miso/

bam_files = [
   "SRR1032173.Aligned.sortedByCoord.out.bam",
   "SRR1032174.Aligned.sortedByCoord.out.bam",
   "SRR1032175.Aligned.sortedByCoord.out.bam",
   "SRR1032176.Aligned.sortedByCoord.out.bam",
   "SRR1032177.Aligned.sortedByCoord.out.bam",
   "SRR1032178.Aligned.sortedByCoord.out.bam"
]

miso_files = [
   "SRR1032173",
   "SRR1032174",
   "SRR1032175",
   "SRR1032176",
   "SRR1032177",
   "SRR1032178"
]

[plotting]
# Dimensions of figure to be plotted (in inches)
fig_width = 7
fig_height = 5
# Factor to scale down introns and exons by
intron_scale = 30
exon_scale = 4
# Whether to use a log scale or not when plotting
logged = False
font_size = 6

# Max y-axis
ymax = 150

# Whether to plot posterior distributions inferred by MISO
show_posteriors = True

# Whether to show posterior distributions as bar summaries
bar_posteriors = False

# Whether to plot the number of reads in each junction
number_junctions = True

resolution = .5
posterior_bins = 40
gene_posterior_ratio = 5

# List of colors for read densities of each sample
colors = [
   "#CC0011",
   "#CC0011",
   "#FF8800",
   "#FF8800",
   "#0080FF",
   "#0080FF",
]
```

```bash
mkdir -p plots/miso

hyperfine --warmup 1 'sashimi_plot --plot-event gene:ENSG00000186092 ref/miso_index/ miso/settings.txt --output-dir plots/miso'

# Time (mean ± σ):      2.400 s ±  0.088 s    [User: 1.608 s, System: 0.702 s]
# Range (min … max):    2.277 s …  2.531 s    10 runs
```


## sashimipy

`sashimipy.txt`
```bash
```

```bash
conda actiavte sashimipy
hyperfine 'sashimipy --density sashimipy.txt -r ref/Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz -e 1:10000-20000 -o test.pdf'

Time (mean ± σ):      4.681 s ±  0.261 s    [User: 9.457 s, System: 10.256 s]
Range (min … max):    4.352 s …  5.047 s    10 runs

Time (mean ± σ):      4.962 s ±  0.238 s    [User: 7.946 s, System: 11.926 s]
Range (min … max):    4.788 s …  5.412 s    10 runs

hyperfine 'sashimipy --density sashimipy.txt -r ref/Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz -e 1:10000-20000 -o test.pdf -p 6'
Time (mean ± σ):      5.026 s ±  0.311 s    [User: 9.898 s, System: 17.662 s]
Range (min … max):    4.566 s …  5.707 s    10 runs
```

## ggsashimi

`ggsashimi.txt`: **Note: must be '\t' seperated**
```bash
HEK293_Control_R1   ./STAR/SRR1032173.Aligned.sortedByCoord.out.bam
HEK293_Control_R2   ./STAR/SRR1032174.Aligned.sortedByCoord.out.bam
Stau1_KnockDown_R1  ./STAR/SRR1032175.Aligned.sortedByCoord.out.bam
Stau1_KnockDown_R2  ./STAR/SRR1032176.Aligned.sortedByCoord.out.bam
Stau1_Overexpression_R1 ./STAR/SRR1032177.Aligned.sortedByCoord.out.bam
Stau1_Overexpression_R2 ./STAR/SRR1032178.Aligned.sortedByCoord.out.bam
```

```bash
conda actiavte ggsashimi
hyperfine 'ggsashimi.py -b ggsashimi.txt -g ref/Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz -c 1:10000-20000 -o test'

Time (mean ± σ):     10.981 s ±  0.458 s    [User: 15.974 s, System: 9.234 s]
Range (min … max):   10.127 s … 11.532 s    10 runs
```

## SplicePlot [failed]

### Prepare SplicePlot

The `SplicePlot-1.1/map_file.txt` as follows:

```bash
HEK293_Control_R1   ../STAR/SRR1032173.Aligned.sortedByCoord.out.bam
HEK293_Control_R2   ../STAR/SRR1032174.Aligned.sortedByCoord.out.bam
Stau1_KnockDown_R1  ../STAR/SRR1032175.Aligned.sortedByCoord.out.bam
Stau1_KnockDown_R2  ../STAR/SRR1032176.Aligned.sortedByCoord.out.bam
Stau1_Overexpression_R1 ../STAR/SRR1032177.Aligned.sortedByCoord.out.bam
Stau1_Overexpression_R2 ../STAR/SRR1032178.Aligned.sortedByCoord.out.bam
```

```bash
conda activate spliceplot

cd SplicePlot-1.1/
zcat ../ref/Homo_sapiens.GRCh38.101.chr.gtf.gz > ../ref/Homo_sapiens.GRCh38.101.chr.gtf
python prepare_annotation.py ../ref/Homo_sapiens.GRCh38.101.chr.gtf ../ref/Homo_sapiens.GRCh38.101.chr.spliceplot.gtf

# we initialized the chr1 for following testings.
python junction_from_window.py 1:1-24895642 map_file.txt > junctions.txt

# we initialized the chr1 for following testings.
# 1:17055-17915,1:17055-17606,1:17055-17233
python initialize_data.py 1:10583 1:18061-18268,1:18061-24738,1:18061-29321,1:18061-195263,1:18061-199837 \
  --vcf ../ref/00-All.vcf.gz \
  --gtf ../ref/Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz \
  --mf map_file.txt
  
mkdir plots
python plot.py pickle_files/1:10583@1:18061-18268,1:18061-24738,1:18061-29321,1:18061-195263,1:18061-199837.p pickle sashimi_test_settings
```

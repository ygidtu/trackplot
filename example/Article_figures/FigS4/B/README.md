## FigS 4B

### An example for visualizing Single cell transcriptional and chromatin accessibility profiling datasets.

Before running this example, you should download these data.

```bash
aria2c -x 8 https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_gex_possorted_bam.bam
aria2c -x 8 https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_gex_possorted_bam.bam.bai
aria2c -x 8 https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam
aria2c -x 8 https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam.bai
aria2c -x 8 https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_cut_sites.bigwig
aria2c -x 8 https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peaks.bed

```

Prepare reference file

```bash
aria2c -c https://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
bedtools sort -i Homo_sapiens.GRCh38.101.chr.gtf.gz | bgzip > Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz
tabix -p gff Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz
```

The cell metadata information,

```bash
#label(must be same to label in density file) cell_barcode  cell_type color_of_cell
ATAC	AAACAGCCAAGGAATC-1	CD4_Naive	#4776a8
ATAC	AAACAGCCAATGCGCT-1	CD4_Naive	#4776a8
ATAC	AAACAGCCAGGATAAC-1	CD4_Naive	#4776a8
ATAC	AAACATGCACTTGTTC-1	CD4_Naive	#4776a8
ATAC	AAACCAACACTAAGAA-1	CD4_Naive	#4776a8
ATAC	AAACCGCGTATGGTGC-1	CD4_Naive	#4776a8
ATAC	AAACCGGCATTAGCCA-1	CD4_Naive	#4776a8
ATAC	AAACGCGCACCTACTT-1	CD4_Naive	#4776a8
```


The command line for generating the plot,

```bash

python ../../../../main.py \
  --barcode cd4_naive_cd16_mono.list \
  --density density.list \
  --dpi 300 --output U2AF1L4.pdf \
  -e chr19:35742064-35746000:+ \
  -r Homo_sapiens.GRCh38.101.sorted.gtf.gz \
  --raster --width 6 -t 1000000 --intron-scale 0.0001 --same-y --sites 35742529,35743011


```

if trackplot was installed by docker, here is the cmd

```bash

cat density.list |grep -v '^#' | while read line; do echo $PWD/${line}; done > density_list_abspath.tsv

docker run -v $PWD:$PWD --rm ygidtu/trackplot  \
  --barcode $PWD/cd4_naive_cd16_mono.list \
  --density $PWD/density_list_abspath.tsv \
  --dpi 300 --output $PWD/U2AF1L4.pdf \
  -e chr19:35742064-35746000:+ \
  -r $PWD/Homo_sapiens.GRCh38.101.sorted.gtf.gz \
  --raster --width 6 -t 1000000 --intron-scale 0.0001 --same-y --sites 35742529,35743011

```

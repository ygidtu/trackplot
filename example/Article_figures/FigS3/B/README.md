## FigS 3B

### An example for visualizing Single cell transcriptional and chromatin accessibility datasets.

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


The command line for generating the plot,

```bash

python ../../../../main.py \
  --barcode cd4_naive_cd16_mono.list \
  --density density.list \
  --dpi 300 --output U2AF1L4.pdf \
  -e chr19:35742064-35746000:+ \
  -r /mnt/raid61/Ref/HomSap/release101/Homo_sapiens.GRCh38.101.sorted.gtf.gz \
  --raster --width 6 -t 1000000 --intron-scale 0.0001 --same-y


```
## FigS 3B

Before running this example, you should download these data.

```bash
aria2c -x 8 https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_gex_possorted_bam.bam
aria2c -x 8 https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_gex_possorted_bam.bam.bai
aria2c -x 8 https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam
aria2c -x 8 https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam.bai
aria2c -x 8 https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_cut_sites.bigwig
aria2c -x 8 https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peaks.bed

aria2c -x 8 https://encode-public.s3.amazonaws.com/2021/10/28/6f0cc163-86c7-4a68-baac-65af90f5a90d/ENCFF053VBX.hic?response-content-disposition=attachment%3B%20filename%3DENCFF053VBX.hic&AWSAccessKeyId=ASIATGZNGCNXVLJBX4VD&Signature=4750a2ET1hFGYPknzGjksVdDQ9M%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEHAaCXVzLXdlc3QtMiJHMEUCIGEw7iFi4aYX%2FW5NabY7DU6Vcj9mpe3eEaFr1ct8MD5KAiEAsYOZi5pJyHRIjB2OHixY0t4y4w6vF9ll8ze2SbrI2y8qzAQIKRAAGgwyMjA3NDg3MTQ4NjMiDPvp9Tj%2FSMgsAy%2B%2FbyqpBFhgsCPcFyyQk1QP2PxwKNgZwFRqWhs4jYmhZ242tLqnoJ%2BoLUPsk556Hul%2FOSiSDIf2JO8CHtvpLCaJXk3sfr%2FUqmgYC1sFVyfL2l0kDMpOznIuhraUEW%2BP0jhrLy%2Bs3L%2F7Fss6TQGaLxWq2AOn6Wgx6m94zbCVniAdMiXus2LsK1lNhzo%2Bs2sZLAsboLfZRiHzwfsVS15QVgL5wtbO1GdFCauPof8u84rATSVYZbTLDmAiF6JESkQ%2F25arZzI%2BheQdDSyG6KhWtpacb2Ey7YOqdTfLoKv9FinvKH3RdVFrOSDBg%2FGiI%2BxvT2ypfVdspnckfpqJv1nGA6niCvLNJRsYr11d8QkYZVFsjFhvEpE6aD%2BtzNzMLCcehVAhwp5soMiuEE9nyG1%2FO3ppbEOk27YLDE4Hty7XMct1vILIHglzbW7dDkW8L6dZe289KKH3JlIo%2FUc5%2B0DxVW19iQ3hJUUofT1aBQ2Eso5U3Elj2cqPNjKUZUDT%2B0DhkFRgxpcjXgsHMiYt2ZbN8mq7VqKyKY0JIwbJmPQRKTMEQcLagWZXYwndd5eisARmjRr4EY%2BEglnmIVnWcUs%2BacR07XD7AkClUbzFVWKLa2MLWJQtrTRVdLc7GpfzfmOo8SI%2FZJk%2F7wzb%2BXypQPFA2JiNOesvQpzqqIkdRSYEFho0%2BHBnCsNArSS8z%2FAC%2FnJbmfz84Gw%2FtO%2BvyTm11WmWJvXI3uxGpzde6PFrbCZKWhswmILQmQY6qQGdsfnaE8ZgvsAyof1TegbVveu8kGSGCuGal9dfJzDWdXeNHwz6Y1zHh%2FGAN%2BuzFNNJ2sNWOsmd9dm39wfQmc3BqmIpcKCDcwunbpxOQZR5G%2F2YrbhZ8myS6Cn16iyWmfvpAWIeXJgDYryBhDcLW%2ByD5kow6tUdDaabdzIS0uHD9ChLwgBQx23laV4FpbXhk0erBTV3a1k%2FIVazpFsh%2F8Tnzn67Tld3s7dq&Expires=1664482722

```

Then you should install [cooler](https://github.com/open2c/cooler), [pairix](https://github.com/4dn-dcic/pairix) and [HiCExplorer](https://github.com/deeptools/HiCExplorer), and convert pairs format into naive Hi-C matrix

```bash
cooler cload pairix -p 16 /mnt/raid61/Ref/HomSap/release101/hg38.chrom.sizes:1000 ENCFF931NQV.pairs.gz ENCFF931NQV_1kb.cool
hicConvertFormat -m ENCFF931NQV_1kb.cool --inputFormat cool --outputFormat h5 -o ENCFF931NQV_1kb.h5

# or just chr19

pairix ENCFF931NQV.pairs.gz 'chr19' | bgzip > ENCFF931NQV.chr19.pairs.gz
cooler cload pairix -p 16 /mnt/raid61/Ref/HomSap/release101/hg38.chrom.sizes:1000 ENCFF931NQV.chr19.pairs.gz ENCFF931NQV.chr19_1kb.cool
hicConvertFormat -m ENCFF931NQV.chr19_1kb.cool --inputFormat cool --outputFormat h5 -o ENCFF931NQV.chr19_1kb.h5

```

Prepare reference file

```bash
aria2c -c https://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
bedtools sort -i Homo_sapiens.GRCh38.101.chr.gtf.gz | bgzip > Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz
tabix -p gff Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz
```


The command line for generating the plot,
```bash

python ../../../..//main.py \
  --hic hic.tsv \
  --barcode Tcell_barcode_sub.list \
  --density density.list \
  --dpi 300 --output scrna.tcell.pdf \
  -e chr19:35600000-35800000:+ \
  -r Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz \
  --raster --width 3 -t 1000000 --choose-primary


```
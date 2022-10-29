
## Figure 2B

Before running this example, you should download these data.

```bash
aria2c -x 16 https://www.encodeproject.org/files/ENCFF125RUG/@@download/ENCFF125RUG.bam
aria2c -x 16 https://www.encodeproject.org/files/ENCFF854PFR/@@download/ENCFF854PFR.bam
aria2c -x 16 https://www.encodeproject.org/files/ENCFF709LHN/@@download/ENCFF709LHN.bam
aria2c -x 16 https://www.encodeproject.org/files/ENCFF613CGT/@@download/ENCFF613CGT.bam

aria2c -x 16 https://www.encodeproject.org/files/ENCFF936SHU/@@download/ENCFF936SHU.bigWig
aria2c -x 16 https://www.encodeproject.org/files/ENCFF363UDO/@@download/ENCFF363UDO.bigWig
aria2c -x 16 https://www.encodeproject.org/files/ENCFF051PIE/@@download/ENCFF051PIE.bigWig
aria2c -x 16 https://www.encodeproject.org/files/ENCFF245YUN/@@download/ENCFF245YUN.bigWig
aria2c -x 16 https://www.encodeproject.org/files/ENCFF556EQK/@@download/ENCFF556EQK.bigWig
aria2c -x 16 https://www.encodeproject.org/files/ENCFF476HFB/@@download/ENCFF476HFB.bigWig

aria2c -x 16 http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons470way/hg38.phastCons470way.bw
```

for generating the plot,

```bash

python sashimi.py/main.py \
  -e chr9:112296343-112335026 \
  -r /mnt/raid61/Ref/HomSap/release101/Homo_sapiens.GRCh38.101.sorted.gtf.gz \
  --density bam.tsv \
  -o PTBP3.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 --show-junction-num \
  --included-junctions chr9:112297918-112330441,chr9:112297918-112333470,chr9:112330476-112333470

```
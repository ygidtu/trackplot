
## FigS 2A

Data was downloaded from [SRX8994511](https://www.ncbi.nlm.nih.gov/sra/SRX8994511). 
After mapping to the corresponding genome using minimap2, we convert bam file into bed file. the detailed command line as follow,

```bash

samtools view -u -F 2304 SRX8994511.example.bam \
| bedtools bamtobed -bed12 -cigar \
| bedtools sort -i /dev/stdin \
| bgzip > SRX8994511.bed.gz

```

Prepare reference file

```bash
aria2c -c https://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
bedtools sort -i Homo_sapiens.GRCh38.101.chr.gtf.gz | bgzip > Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz
tabix -p gff Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz
```


the command line for generating the plots.

```bash

python  ../../../../main.py \
  -r Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz \
  -e 21:43092956-43107570:+ \
  --density bam.tsv \
  --igv igv.tsv \
  --focus 43100453-43100519:43101366-43101432 \
  -o igv_plot.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 --show-junction-num --domain

```
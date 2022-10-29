
## Figure 3A

Data was downloaded from [SRX8994511](https://www.ncbi.nlm.nih.gov/sra/SRX8994511). 
After mapping to the corresponding genome using minimap2, we convert bam file into bed file. the detailed command line as follow,

```bash

samtools view -u -F 2304 SRX8994511.example.bam \
| bedtools bamtobed -bed12 -cigar \
|bedtools sort -i /dev/stdin \
| bgzip > SRX8994511.bed.gz

```

the command line for generating the plots.

```bash

python  sashimi.py/main.py \
  -r /mnt/raid61/Ref/HomSap/release101/Homo_sapiens.GRCh38.101.sorted.gtf.gz \
  -e 21:43092956-43107570:+ \
  --density bam.tsv \
  --igv igv.tsv \
  --focus 43100453-43100519:43101366-43101432 \
  -o igv_plot.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 --show-junction-num --domain 

```

## FigS 1A 

Data was downloaded from [Wang et al, Cell Research, 2021.](https://www.nature.com/articles/s41422-020-00451-z)

Prepare reference file

```bash
aria2c -c https://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
bedtools sort -i Homo_sapiens.GRCh38.101.chr.gtf.gz | bgzip > Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz
tabix -p gff Homo_sapiens.GRCh38.101.chr.sorted.gtf.gz
```


Run this example,

```bash

python ../../../../main.py \
  -e 7:107830089-107832297 \
  -r ./Homo_sapiens.GRCh37.87.gtf.sorted.gz \
  --density bam.tsv \
  -o TNP.NRCAM.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 --show-junction-num -t 10 --domain

```

if trackplot was installed by docker, here is the cmd

```bash

cat bam.tsv | while read line; do echo $PWD/${line}; done > bam_abspath.tsv

docker run -v $PWD:$PWD --rm ygidtu/trackplot -e 7:107830089-107832297 \
  -r $PWD/Homo_sapiens.GRCh37.87.gtf.sorted.gz \
  --density $PWD/bam_abspath.tsv \
  -o $PWD/TNP.NRCAM.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 --show-junction-num -t 10 --domain

```

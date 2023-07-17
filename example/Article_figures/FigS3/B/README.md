
## FigS 3B

Data was download from [Gao et al. Genome Biology,2021.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02241-7),
and we used [nanom6A](https://github.com/gaoyubang/nanom6A) and [nanopolish](https://github.com/nanoporetech/pipeline-polya-ng) for calling m6A and length of poly(A)  

Prepare reference file

```bash
aria2c -c https://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz
bedtools sort -i Homo_sapiens.GRCh37.87.chr.gtf.gz | bgzip > Homo_sapiens.GRCh37.87.chr.sorted.gtf.gz
tabix -p gff Homo_sapiens.GRCh37.87.chr.sorted.gtf.gz
```


```bash
python ../../../../main.py \
    -r Homo_sapiens.GRCh37.87.gtf.sorted.gz \
    -e chr7:5566600-5570232 \
    --igv igv.tsv \
    -o igv.ATCB.pdf \
    --dpi 300 \
    --width 10 \
    --height 0.5 --rs rs --polya pa --m6a ma

```

if trackplot was installed by docker, here is the cmd

```bash

cat igv.tsv |grep -v '^#' | while read line; do echo $PWD/${line}; done > igv_abspath.tsv

docker run -v $PWD:$PWD --rm ygidtu/trackplot  \
    -r $PWD/Homo_sapiens.GRCh37.87.gtf.sorted.gz \
    -e chr7:5566600-5570232 \
    --igv $PWD/igv_abspath.tsv \
    -o $PWD/igv.ATCB.pdf \
    --dpi 300 \
    --width 10 \
    --height 0.5 --rs rs --polya pa --m6a ma

```


## FigS 2B

Data was download from [Gao et al. Genome Biology,2021.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02241-7),
and we used [nanom6A](https://github.com/gaoyubang/nanom6A) and [nanopolish](https://github.com/nanoporetech/pipeline-polya-ng) for calling m6A and length of poly(A)  


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
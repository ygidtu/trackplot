
## FigS 3A

Data was downloaded from [Zhou et al., Nucleic Acids Research, 2022](https://academic.oup.com/nar/article/50/11/e66/6548409)

```bash
aria2c -c https://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.chr.gtf.gz
bedtools sort -i Mus_musculus.GRCm38.101.chr.gtf.gz | bgzip > Mus_musculus.GRCm38.101.chr.sorted.gtf.gz
tabix -p gff Mus_musculus.GRCm38.101.chr.sorted.gtf.gz
```

The command line for generating the plot
```bash

python ../../../..//main.py \
  -e 17:35832921-35835600 \
  -r Mus_musculus.GRCm38.101.chr.sorted.gtf.gz \
  --density density_list.tsv \
  -o hsc_8w.Tubb5.remove_dup.2.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 -t 100000 \
  --barcode cell_meta.tsv --remove-duplicate-umi -p 12


```

if trackplot was installed by docker, here is the cmd

```bash

cat density_list.tsv |grep -v '^#' | while read line; do echo $PWD/${line}; done > density_list_abspath.tsv

docker run -v $PWD:$PWD --rm ygidtu/trackplot  \
  -e 17:35832921-35835600 \
  -r $PWD/Mus_musculus.GRCm38.101.chr.sorted.gtf.gz \
  --density $PWD/density_list_abspath.tsv \
  -o $PWD/hsc_8w.Tubb5.remove_dup.2.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 -t 100000 \
  --barcode $PWD/cell_meta.tsv --remove-duplicate-umi -p 12

```

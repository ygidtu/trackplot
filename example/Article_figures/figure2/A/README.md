
## Figure 2A 

Data was downloaded from [Wang et al, Cell Research, 2021.](https://www.nature.com/articles/s41422-020-00451-z)

Run this example,

```bash

python sashimi.py/main.py \
  -e 7:107830089-107832297 \
  -r /mnt/raid61/public_data/transferdata/GBM_bulk/data/ref/Homo_sapiens.GRCh37.87.gtf.sorted.gz \
  --density bam.tsv \
  -o TNP.NRCAM.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 --show-junction-num -t 10 --domain

```
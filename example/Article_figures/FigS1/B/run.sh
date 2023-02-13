python sashimi/sashimi.py/main.py \
  -e chr9:112296343-112335026 \
  -r Homo_sapiens.GRCh38.101.sorted.gtf.gz \
  --density bam.tsv \
  -o PTBP3.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 --show-junction-num \
  --included-junctions chr9:112297918-112330441,chr9:112297918-112333470,chr9:112330476-112333470


python sashimi.py/main.py \
  --hic hic.tsv \
  --barcode Tcell_barcode_sub.list \
  --density density.list \
  --dpi 300 --output scrna.tcell.pdf \
  -e chr19:35600000-35800000:+ \
  -r Homo_sapiens.GRCh38.101.sorted.gtf.gz \
  --raster --width 3 -t 1000000 --choose-primary


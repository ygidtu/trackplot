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
python sashimi.py/main.py \
  -e 17:35832921-35835613 \
  -r Mus_musculus.GRCm38.101.gtf.gz \
  --density density_list.tsv \
  -o hsc_8w.Tubb5.remove_dup.2.pdf \
  --dpi 300 \
  --width 6 \
  --height 1 -t 100000 \
  --barcode cell_meta.tsv --remove-duplicate-umi -p 12


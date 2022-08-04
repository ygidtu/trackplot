# sashimi

## Example

```bash
python main.py \
  -e chr1:1270656-1284730:+ \
  -r example/example.sorted.gtf.gz \
  --interval example/interval_list.tsv \
  --density example/density_list.tsv \
  --show-side \
  --igv example/igv.tsv \
  --heatmap example/heatmap_list.tsv \
  --focus 1272656-1272656:1275656-1277656 \
  --stroke 1275656-1277656:1277856-1278656@blue \
  --sites 1271656,1271656,1272656 \
  --line example/line_list.tsv \
  -o example/example.png \
  --dpi 300 \
  --width 10 \
  --height 1 \
  --barcode example/barcode_list.tsv \
  --domain
```

![](example/example.png)


mkdir -p out

../fast_bn --input data_all/all_disc.tsv --score bdeu \
  --ess 10 --tabu 30 --iters 5000 --topk 20 \
  --jindex-cache 1024 \
  --save        out/edges.tsv \
  --save-names  out/edges_named.tsv \
  --save-counts out/all_counts.tsv >out/log.txt 2>&1


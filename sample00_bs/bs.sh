../fast_bn --input ./all_disc10.tsv --score bic \
  --bootstrap 100 \
  --seed 42 \
  --save-bootstrap-counts boot_edges.tsv \
  --tabu 20 --iters 3000 > log.txt 2>&1

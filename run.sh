#./fast_bn all_disc10.tsv bdeu --ess 10 --tabu 30 --iters 8000
#./fast_bn all_disc1000.tsv bic --tabu 30 --iters 8000
./fast_bn --input  all_disc10.tsv --score bdeu --ess 10 --tabu 30 --iters 8000 \
  --save init_edges.tsv \
  --save-names init_edges_named.tsv \
  --save-counts all_counts.tsv


#./fast_bn all_disc10.tsv bdeu --ess 10 --tabu 30 --iters 8000
#./fast_bn all_disc1000.tsv bic --tabu 30 --iters 8000
#./fast_bn all_disc.tsv bdeu --ess 10 --tabu 30 --iters 1000 \
#  --save init_edges.tsv \
#  --save-names init_edges_named.tsv \
#  --save-counts all_counts.tsv
./fast_bn dummy bic \
  --score-dataset all_disc.tsv \
  --init init_edges.tsv \
  --counts all_counts.tsv \
  --alpha 1.0


./fast_bn dummy bic \
  --edge-importance \
  --score-dataset all_disc.tsv \
  --init init_edges.tsv \
  --counts all_counts.tsv \
  --alpha 1.0 \
  --ess 10.0 \
  --save-edge-importance edge_importance.tsv


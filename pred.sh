#./fast_bn all_disc10.tsv bdeu --ess 10 --tabu 30 --iters 8000
#./fast_bn all_disc1000.tsv bic --tabu 30 --iters 8000
#./fast_bn all_disc.tsv bdeu --ess 10 --tabu 30 --iters 1000 \
#  --save init_edges.tsv \
#  --save-names init_edges_named.tsv \
#  --save-counts all_counts.tsv
./fast_bn --input dummy --score bic \
  --score-dataset data/all_disc10.tsv \
  --init out/init_edges.tsv \
  --counts out/all_counts.tsv \
  --alpha 1.0


./fast_bn --input dummy --score bic \
  --edge-importance \
  --score-dataset data/all_disc10.tsv \
  --init out/init_edges.tsv \
  --counts out/all_counts.tsv \
  --alpha 1.0 \
  --ess 10.0 \
  --save-edge-importance out/edge_importance.tsv


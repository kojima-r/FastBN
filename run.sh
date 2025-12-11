#./fast_bn all_disc10.tsv bdeu --ess 10 --tabu 30 --iters 8000
#./fast_bn all_disc1000.tsv bic --tabu 30 --iters 8000
mkdir -p out
./fast_bn --input  data/all_disc10.tsv --score bdeu --ess 10 --tabu 30 --iters 8000 \
  --save out/init_edges.tsv \
  --save-names out/init_edges_named.tsv \
  --save-counts out/all_counts.tsv


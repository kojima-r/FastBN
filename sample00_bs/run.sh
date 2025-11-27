./fast_bn ./all_disc.tsv bdeu --bootstrap 100  --bootstrap-seed 1 --save-bootstrap-counts boot_edges.tsv  --tabu 20 --iters 3000  >log1.txt &
./fast_bn ./all_disc.tsv bdeu --bootstrap 100  --bootstrap-seed 2 --save-bootstrap-counts boot_edges.tsv  --tabu 20 --iters 3000  >log2.txt &
./fast_bn ./all_disc.tsv bdeu --bootstrap 100  --bootstrap-seed 3 --save-bootstrap-counts boot_edges.tsv  --tabu 20 --iters 3000  >log3.txt &
./fast_bn ./all_disc.tsv bdeu --bootstrap 100  --bootstrap-seed 4 --save-bootstrap-counts boot_edges.tsv  --tabu 20 --iters 3000  >log4.txt &
./fast_bn ./all_disc.tsv bdeu --bootstrap 100  --bootstrap-seed 5 --save-bootstrap-counts boot_edges.tsv  --tabu 20 --iters 3000  >log5.txt &

wait


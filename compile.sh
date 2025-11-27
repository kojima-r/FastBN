g++ -O3 -march=native -std=c++17 fast_bn.cpp -o fast_bn
#time ./fast_bn --input all_disc100.tsv --score bdeu --ess 10 --tabu 30 --iters 100


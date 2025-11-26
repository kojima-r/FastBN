g++ -pg -std=c++17 fast_bn.cpp -o fast_bn_pg

#./fast_bn_pg --input all_disc100.tsv --score bdeu --ess 10 --tabu 30 --iters 100
#gprof ./fast_bn_pg gmon.out

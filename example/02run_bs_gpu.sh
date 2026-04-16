#!/bin/bash

export NV_ACC_TIME=1

bin=../fast_bn
output=./`date +%Y%m%d_%H%M%S`_02_tri
mkdir -p ${output}

for seed in `seq 1 5`
do

${bin} --init ./out/edges.tsv \
    --input ./data_bin/all_disc100.tsv \
    --score bdeu --bootstrap 10 \
    --save-bootstrap-counts bs/edges.tsv \
    --topk 20 --jindex-cache 1024 --tabu 20 --iters 5000 \
    --seed ${seed} > ${output}/02log_bs${seed}.txt 2>&1

python check_result.py ${output}/edges_bs_seed000${seed}.tsv out/edges_bs_seed000${seed}.tsv

done

cp ${bin} ${output}

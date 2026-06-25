#!/bin/bash

export NV_ACC_TIME=1
#export NVCOMPILER_ACC_DEBUG=1

bin=../fast_bn
output=./`date +%Y%m%d_%H%M%S`_gpu_01_tri
mkdir -p ${output}

${bin} --input ./data_tri/all_disc_tri100.tsv --score bdeu \
  --ess 10 --tabu 30 --iters 5000 --topk 20 \
  --jindex-cache 1024 \
  --save        ${output}/edges_tri.tsv \
  --save-names  ${output}/edges_named_tri.tsv \
  --save-counts ${output}/all_counts_tri.tsv > ${output}/01log_tri.txt 2>&1

python check_result.py ${output}/edges_tri.tsv out/edges_tri.tsv
python check_result.py ${output}/edges_named_tri.tsv out/edges_named_tri.tsv
python check_result.py ${output}/all_counts_tri.tsv out/all_counts_tri.tsv

cp ${bin} ${output}

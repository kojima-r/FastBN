#!/bin/bash

export NV_ACC_TIME=1
#export NVCOMPILER_ACC_DEBUG=1

bin=../fast_bn
output=./`date +%Y%m%d_%H%M%S`_gpu_03_bs
mkdir -p ${output}

# Run without searching (iters=0) to output integ_all_counts.tsv
${bin} --input ./data_bin/all_disc100.tsv --score bdeu \
  --init out/integ_edges.tsv\
  --ess 10 --tabu 30 --iters 0 --topk 20 \
  --jindex-cache 1024 \
  --save        ${output}/integ_edges2.tsv \
  --save-names  ${output}/integ_edges_named.tsv \
  --save-counts ${output}/integ_all_counts.tsv > ${output}/03log_bs.txt 2>&1

python check_result.py ${output}/integ_edges2.tsv out/integ_edges2.tsv
python check_result.py ${output}/integ_edges_named.tsv out/integ_edges_named.tsv
python check_result.py ${output}/integ_all_counts.tsv out/integ_all_counts.tsv

cp ${bin} ${output}

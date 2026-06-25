#!/bin/bash

export NV_ACC_TIME=1
#export NVCOMPILER_ACC_DEBUG=1

bin=../fast_bn
output=./`date +%Y%m%d_%H%M%S`_gpu_04
mkdir -p ${output}

${bin} --score bic \
  --score-dataset ./data_bin/all_disc100.tsv \
  --init out/integ_edges_score.tsv \
  --counts out/integ_all_counts.tsv \
  --alpha 1.0 > ${output}/04score_check.txt 2>&1

cp ${bin} ${output}

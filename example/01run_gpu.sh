#!/bin/bash
#
# 単純な探索モード
#

export NV_ACC_TIME=1
export NVCOMPILER_ACC_NOTIFY=3
#export NVCOMPILER_ACC_DEBUG=0x800

bin=../fast_bn
input=./data_bin/all_disc100.tsv
output=./`date +%Y%m%d_%H%M%S`_01
mkdir -p ${output}

# profile を取る時は NV_ACC_TIME を無効にすること
# compute-sanitizer \
# nsys profile -t cuda,openacc -o {output}/profile_output \
${bin} --input ${input} --score bdeu \
  --ess 10 --tabu 30 --iters 5000 --topk 20 \
  --jindex-cache 1024 \
  --save        ${output}/edges.tsv \
  --save-names  ${output}/edges_named.tsv \
  --save-counts ${output}/all_counts.tsv > ${output}/01log.txt 2>&1

python3 check_result.py ${output}/edges.tsv out/edges.tsv
python3 check_result.py ${output}/edges_named.tsv out/edges_named.tsv
python3 check_result.py ${output}/all_counts.tsv out/all_counts.tsv

cp ${bin} ${output}

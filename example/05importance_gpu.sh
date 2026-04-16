#!/bin/bash

export NV_ACC_TIME=1

bin=../fast_bn
output=./`date +%Y%m%d_%H%M%S`_05_importance

${bin} --score bic \
  --edge-importance \
  --score-dataset ./data_bin/all_disc100.tsv \
  --init out/integ_edges_score.tsv \
  --counts out/integ_all_counts.tsv \
  --alpha 1.0 \
  --ess 10.0 \
  --save-edge-importance ${output}/edge_importance.tsv

python check_result.py ${output}/edge_importance.tsv out/edge_importance.tsv -t nearly

cp ${bin} ${output}

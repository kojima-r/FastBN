bin=../fast_bn
input=./data_bin/all_disc100.tsv
output=./out

${bin} --score bic \
  --edge-importance \
  --score-dataset ${input} \
  --init out/integ_edges_score.tsv \
  --counts out/integ_all_counts.tsv \
  --alpha 1.0 \
  --ess 10.0 \
  --save-edge-importance ${output}/edge_importance.tsv


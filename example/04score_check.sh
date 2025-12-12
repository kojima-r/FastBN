bin=../fast_bn
input=./data_bin/all_disc100.tsv
output=./out

${bin} --score bic \
  --score-dataset ${input} \
  --init out/integ_edges_score.tsv \
  --counts out/integ_all_counts.tsv \
  --alpha 1.0


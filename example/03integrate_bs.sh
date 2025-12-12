#python ../compute_bs_prob.py  \
#	--input bs/*.tsv \
#	--out-edge out/integ_edges.tsv \
#	--out out/integ_edges_score.tsv

bin=../fast_bn
input=./data_bin/all_disc100.tsv
output=./out
mkdir -p ${output}
# Run without searching (iters=0) to output integ_all_counts.tsv
${bin} --input ${input} --score bdeu \
  --init out/integ_edges.tsv\
  --ess 10 --tabu 30 --iters 0 --topk 20 \
  --jindex-cache 1024 \
  --save        ${output}/integ_edges2.tsv \
  --save-names  ${output}/integ_edges_named.tsv \
  --save-counts ${output}/integ_all_counts.tsv #>out/log.txt 2>&1


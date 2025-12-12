bin=../fast_bn
input=./data_tri/all_disc_tri100.tsv
output=./out_tri
mkdir -p ${output}

${bin} --input ${input} --score bdeu \
  --ess 10 --tabu 30 --iters 5000 --topk 20 \
  --jindex-cache 1024 \
  --save        ${output}/edges.tsv \
  --save-names  ${output}/edges_named.tsv \
  --save-counts ${output}/all_counts.tsv #>out/log.txt 2>&1


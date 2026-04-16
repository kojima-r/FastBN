#!/bin/bash

export NV_ACC_TIME=1

bin=../fast_bn
output=./`date +%Y%m%d_%H%M%S`
mkdir -p ${output}

${bin} --input ./data_bin/all_disc100.tsv --score bdeu \
  --ess 10 --tabu 30 --iters 5000 --topk 20 \
  --jindex-cache 1024 \
  --save        ${output}/edges.tsv \
  --save-names  ${output}/edges_named.tsv \
  --save-counts ${output}/all_counts.tsv > ${output}/01log.txt 2>&1

python check_result.py ${output}/edges.tsv out/edges.tsv
python check_result.py ${output}/edges_named.tsv out/edges_named.tsv
python check_result.py ${output}/all_counts.tsv out/all_counts.tsv

${bin} --input ./data_tri/all_disc_tri100.tsv --score bdeu \
  --ess 10 --tabu 30 --iters 5000 --topk 20 \
  --jindex-cache 1024 \
  --save        ${output}/edges_tri.tsv \
  --save-names  ${output}/edges_named_tri.tsv \
  --save-counts ${output}/all_counts_tri.tsv > ${output}/01log_tri.txt 2>&1

python check_result.py ${output}/edges_tri.tsv out/edges_tri.tsv
python check_result.py ${output}/edges_named_tri.tsv out/edges_named_tri.tsv
python check_result.py ${output}/all_counts_tri.tsv out/all_counts_tri.tsv

for seed in `seq 1 5`
do

${bin} --init ./out/edges.tsv \
    --input ./data_bin/all_disc100.tsv \
    --score bdeu --bootstrap 10 \
    --save-bootstrap-counts ${output}/edges_bs.tsv \
    --topk 20 --jindex-cache 1024 --tabu 20 --iters 5000 \
    --seed ${seed} > ${output}/02log_bs${seed}.txt 2>&1

python check_result.py ${output}/edges_bs_seed000${seed}.tsv out/edges_bs_seed000${seed}.tsv

done

for seed in `seq 1 5`
do

${bin} --init ./out/edges_tri.tsv \
    --input ./data_tri/all_disc_tri100.tsv \
    --score bdeu --bootstrap 10 \
    --save-bootstrap-counts ${output}/edges_bs_tri.tsv \
    --topk 20 --jindex-cache 1024 --tabu 20 --iters 5000 \
    --seed ${seed} > ${output}/02log_bs_tri${seed}.txt 2>&1

python check_result.py ${output}/edges_bs_tri_seed000${seed}.tsv out/edges_bs_tri_seed000${seed}.tsv

done

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

${bin} --score bic \
  --score-dataset ./data_bin/all_disc100.tsv \
  --init out/integ_edges_score.tsv \
  --counts out/integ_all_counts.tsv \
  --alpha 1.0 > ${output}/04score_check.txt 2>&1

${bin} --score bic \
  --edge-importance \
  --score-dataset ./data_bin/all_disc100.tsv \
  --init out/integ_edges_score.tsv \
  --counts out/integ_all_counts.tsv \
  --alpha 1.0 \
  --ess 10.0 \
  --save-edge-importance ${output}/edge_importance.tsv> ${output}/05importance_check.txt 2>&1

python check_result.py ${output}/edge_importance.tsv out/edge_importance.tsv -t nearly

cp ${bin} ${output}

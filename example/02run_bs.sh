bin=../fast_bn
input=./data_bin/all_disc100.tsv
output=./bs
mkdir -p ${output}

for seed in `seq 1 5`
do

${bin} --init ./out/edges.tsv \
	--input ${input} \
	--score bdeu --bootstrap 10 \
       	--save-bootstrap-counts bs/edges.tsv \
	--topk 20 --jindex-cache 1024 --tabu 20 --iters 5000 \
	--seed ${seed} >bs/log${seed}.txt 2>&1 &

done

wait


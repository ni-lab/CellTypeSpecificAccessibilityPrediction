#!/bin/bash
imbalance_dir=$1
genome=$2
targets=$3
params=$4
model=$5
model_dir=${model::-9}
output=$6
output_dir=${output::-27}

echo ${model_dir}
echo ${output_dir}


mkdir ${output_dir}

for cell_type in `ls -1 ${imbalance_dir}`;
do
    mkdir ${output_dir}/${cell_type}

    for cluster in `ls -1 ${imbalance_dir}/${cell_type} | cut -d '.' -f1`;
    do
    	python /global/home/users/poojakathail/basenji/bin/basenji_sad.py --stats SAD,REF,ALT -f /clusterfs/nilah/pooja/genomes/${genome}.ml.fa -o ${output_dir}/${cell_type}/${cluster} --rc --shifts "1,0,-1" -t ${targets} ${params} ${model_dir}/model_best.h5 ${imbalance_dir}/${cell_type}/${cluster}.vcf
    done
done  

touch ${output}

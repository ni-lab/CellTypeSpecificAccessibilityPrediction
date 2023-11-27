#!/bin/bash

cluster_dir=$1
cluster_list=$2
genome=$3
targets=$4
params=$5
model=$6
model_dir=${model::-9}
output=$7
output_dir=${output::-22}

echo ${model_dir}
echo ${output_dir}

Field_Separator=$IFS
IFS=,

mkdir ${output_dir}
for cluster in ${cluster_list};
do
	python /global/home/users/poojakathail/basenji/bin/basenji_predict_bed.py -f /clusterfs/nilah/pooja/genomes/${genome}.ml.fa -g /clusterfs/nilah/pooja/genomes/human.${genome}.genome -o ${output_dir}/${cluster} --rc --shifts 1,0,-1 -t ${targets} ${params} ${model_dir}/model_best.h5 ${cluster_dir}/${cluster}_test_chrs.bed
done
IFS=$Field_Separator

touch ${output}


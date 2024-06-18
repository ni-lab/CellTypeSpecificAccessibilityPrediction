#!/bin/bash

bed_file_dir=$1
out_dir=$2
cluster_list=$3
g_dir="/clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink"
snp_list="/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/list.txt"

for cluster in ${cluster_list};
do
    parallel -j 11 python /clusterfs/nilah/pooja/software/ldsc/make_annot.py \
        --bed-file ${bed_file_dir}/${cluster}.hg19.bed \
        --bimfile ${g_dir}/1000G.EUR.QC.{}.bim \
        --annot-file ${out_dir}/${cluster}_baseline_snp_list.{}.annot.gz ::: {1..22}
        
    parallel -j 11 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --l2 --bfile ${g_dir}/1000G.EUR.QC.{} --ld-wind-cm 1 --thin-annot --annot ${out_dir}/${cluster}_baseline_snp_list.{}.annot.gz --out ${out_dir}/${cluster}_baseline_snp_list.{} --print-snps ${snp_list} ::: {1..22}
done


#!/bin/bash

dataset=$1
tasks=$2

annot_dir="/clusterfs/nilah/pooja/${dataset}/ldsc/annots"
if [ "${dataset}" = "kidney_data" ]
then
    tissue="Kidney"
else
    tissue="Blood_Immune"
fi

for col in ${tissue}_more_cell_type_specific ${tissue}_less_cell_type_specific;
do
    for sad_stat in max mean;
    do
        for bin in bottom top;
        do
            percentile_bin="0.5"
            annot_col="${col}_${dataset}_${tasks}_${sad_stat}_sad_bin_${bin}_${percentile_bin}"

            parallel -j 22 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --l2 --bfile /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink/1000G.EUR.QC.{} --ld-wind-cm 1 --annot ${annot_dir}/${annot_col}.{}.annot --thin-annot --out ${annot_dir}/${annot_col}.{} --print-snps /clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/list.txt ::: {1..22}
        done
    done
done


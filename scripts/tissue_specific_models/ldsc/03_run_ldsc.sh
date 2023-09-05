#!/bin/bash

dataset=$1
ld_scores_dir=$2
out_dir=$3

if [ $dataset -eq "kidney"]
then
	sum_stats_dir="/clusterfs/nilah/pooja/kidney_data/ldsc/sumstats"
	gwas="kidney_30700_UKB"
else
	sum_stats_dir="/clusterfs/nilah/pooja/immune_atlas/ldsc/sumstats"
	gwas="PASS_Crohns_Disease PASS_Rheumatoid_Arthritis PASS_Ulcerative_Colitis"
fi

ref_ld_list="/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/baselineLD."
for cluster in cluster_cell_type_specific cluster_open;
do
    ref_ld_list="${ref_ld_list},${ld_scores_dir}/${cluster}_baseline_snp_list."
done

parallel -j 3 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --h2 ${sum_stats_dir}/{}.sumstats --ref-ld-chr ${ref_ld_list} --w-ld-chr /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --frqfile-chr /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_frq/1000G.EUR.QC. --overlap-annot --out ${out_dir}/{}_baseline_cts_open_peak_sets ::: ${gwas}


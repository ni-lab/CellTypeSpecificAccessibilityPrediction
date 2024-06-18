#!/bin/bash

dataset=$1
ld_scores_dir=$2
out_dir=$3

if [[ "$dataset" == "kidney" ]]
then
    sum_stats_dir="/clusterfs/nilah/pooja/kidney_data/ldsc/sumstats"
    gwas="kidney_30700_UKBB"
    clusters="all_peaks_non_ubiq cluster3_Ubiquitous cluster11_PT_baseline_snp_list cluster4_DistalNephron_baseline_snp_list cluster9_Stroma_baseline_snp_list"
else
    sum_stats_dir="/clusterfs/nilah/pooja/immune_atlas/ldsc/sumstats/lupus_paper"
    gwas="PASS_Crohns_Disease PASS_Rheumatoid_Arthritis PASS_Ulcerative_Colitis"
    clusters="cluster_cell_type_specific_baseline_snp_list cluster_open_baseline_snp_list cluster_myeloid_resting_baseline_snp_list cluster_T_resting_baseline_snp_list cluster_nk_resting_baseline_snp_list"
fi

for cluster in ${clusters}
do
    ref_ld_list="/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/baselineLD_remove_enhancer_promoter_annots."
    ref_ld_list="${ref_ld_list},${ld_scores_dir}/${cluster}."
    parallel -j 3 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --h2 ${sum_stats_dir}/{}.sumstats --ref-ld-chr ${ref_ld_list} --w-ld-chr /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --frqfile-chr /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_frq/1000G.EUR.QC. --overlap-annot --out ${out_dir}/{}_baseline_remove_enhancer_promoter_annots_${cluster} ::: ${gwas}
done




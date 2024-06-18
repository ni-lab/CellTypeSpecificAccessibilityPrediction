#!/bin/bash

enformer_dir="/global/scratch/users/poojakathail/enformer"
baseline_annot="/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/baselineLD_remove_enhancer_promoter_annots."
weights="/clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
freqs="/clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_frq/1000G.EUR.QC."

traits=("50_irnt" "23104_irnt" "6152_8" "2443" "6152_9" "20116_0" "78_irnt")
tissues=("Musculoskeletal_Connective" "CNS" "Blood_Immune" "Pancreas" "Blood_Immune" "CNS" "Cardiovascular")

dataset="kidney_data"
tissue="Kidney"
ldsc_dir="/clusterfs/nilah/pooja/${dataset}/ldsc"

trait="30700_irnt"
for peaks in less_cell_type_specific more_cell_type_specific;
do
    parallel -j 6 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --h2 ${enformer_dir}/ldsc/sumstats/${trait}.ldsc.imputed_v3.both_sexes.tsv.bgz --ref-ld-chr ${baseline_annot},${ldsc_dir}/annots/${tissue}_${peaks}_${dataset}_{}. --w-ld-chr ${weights} --frqfile-chr ${freqs} --overlap-annot --print-coefficients --out ${ldsc_dir}/partitioned_heritability/${trait}_baseline_remove_enhancer_promoter_annots_${tissue}_${peaks}_${dataset}_{} ::: multitask_max_sad_bin_bottom_0.5 multitask_max_sad_bin_top_0.5 multitask_mean_sad_bin_bottom_0.5 multitask_mean_sad_bin_top_0.5 single_task_max_sad_bin_bottom_0.5 single_task_max_sad_bin_top_0.5 single_task_mean_sad_bin_bottom_0.5 single_task_mean_sad_bin_top_0.5 multitask_8x_params_same_layers_max_sad_bin_bottom_0.5 multitask_8x_params_same_layers_max_sad_bin_top_0.5 multitask_8x_params_same_layers_mean_sad_bin_bottom_0.5 multitask_8x_params_same_layers_mean_sad_bin_top_0.5 
done

dataset="immune_atlas"
tissue="Blood_Immune"
ldsc_dir="/clusterfs/nilah/pooja/${dataset}/ldsc"

for trait in 6152_8 6152_9;
do
    for peaks in less_cell_type_specific more_cell_type_specific;
    do
        parallel -j 6 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --h2 ${enformer_dir}/ldsc/sumstats/${trait}.ldsc.imputed_v3.both_sexes.tsv.bgz --ref-ld-chr ${baseline_annot},${ldsc_dir}/annots/${tissue}_${peaks}_${dataset}_{}. --w-ld-chr ${weights} --frqfile-chr ${freqs} --overlap-annot --print-coefficients --out ${ldsc_dir}/partitioned_heritability/${trait}_baseline_remove_enhancer_promoter_annots_${tissue}_${peaks}_${dataset}_{} ::: multitask_max_sad_bin_bottom_0.5 multitask_max_sad_bin_top_0.5 multitask_mean_sad_bin_bottom_0.5 multitask_mean_sad_bin_top_0.5 single_task_max_sad_bin_bottom_0.5 single_task_max_sad_bin_top_0.5 single_task_mean_sad_bin_bottom_0.5 single_task_mean_sad_bin_top_0.5 multitask_8x_params_same_layers_max_sad_bin_bottom_0.5 multitask_8x_params_same_layers_max_sad_bin_top_0.5 multitask_8x_params_same_layers_mean_sad_bin_bottom_0.5 multitask_8x_params_same_layers_mean_sad_bin_top_0.5 
    done
done
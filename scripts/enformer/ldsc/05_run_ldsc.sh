#!/bin/bash

enformer_dir="/global/scratch/users/poojakathail/enformer"
dummy_annot="/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/baselineLD_dummy_var."
baseline_annot="/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/baselineLD."
baseline_remove_enhancer_promoter_annot="/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/baselineLD_remove_enhancer_promoter_annots."
weights="/clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
freqs="/clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_frq/1000G.EUR.QC."

traits=("50_irnt" "23104_irnt" "6152_8" "2443" "6152_9" "20116_0" "78_irnt", "30700_irnt")
tissues=("Musculoskeletal_Connective" "CNS" "Blood_Immune" "Pancreas" "Blood_Immune" "CNS" "Cardiovascular", "Kidney")

for i in {0..7}
do
    echo ${traits[$i]}
    echo ${tissues[$i]}
    
    # no baseline
    parallel -j 8 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --h2 ${enformer_dir}/ldsc/sumstats/${traits[$i]}.ldsc.imputed_v3.both_sexes.tsv.bgz --ref-ld-chr ${dummy_annot},${enformer_dir}/ldsc/annots/${tissues[$i]}_{}. --w-ld-chr ${weights} --frqfile-chr ${freqs} --overlap-annot --print-coefficients --out ${enformer_dir}/ldsc/out/${traits[$i]}_no_baseline_${tissues[$i]}_{} ::: less_cell_type_specific less_cell_type_specific_mean_sad_bin_bottom_0.5 less_cell_type_specific_mean_sad_bin_top_0.5 more_cell_type_specific more_cell_type_specific_mean_sad_bin_bottom_0.5 more_cell_type_specific_mean_sad_bin_top_0.5 less_cell_type_specific_mean_sad_bin_bottom_0.3 less_cell_type_specific_mean_sad_bin_top_0.3 less_cell_type_specific_mean_sad_bin_bottom_0.1 less_cell_type_specific_mean_sad_bin_top_0.1 more_cell_type_specific_mean_sad_bin_bottom_0.3 more_cell_type_specific_mean_sad_bin_top_0.3 more_cell_type_specific_mean_sad_bin_bottom_0.1 more_cell_type_specific_mean_sad_bin_top_0.1
     
    # baseline
    parallel -j 8 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --h2 ${enformer_dir}/ldsc/sumstats/${traits[$i]}.ldsc.imputed_v3.both_sexes.tsv.bgz --ref-ld-chr ${baseline_annot},${enformer_dir}/ldsc/annots/${tissues[$i]}_{}. --w-ld-chr ${weights} --frqfile-chr ${freqs} --overlap-annot --print-coefficients --out ${enformer_dir}/ldsc/out/${traits[$i]}_baseline_${tissues[$i]}_{} ::: less_cell_type_specific less_cell_type_specific_mean_sad_bin_bottom_0.5 less_cell_type_specific_mean_sad_bin_top_0.5 more_cell_type_specific more_cell_type_specific_mean_sad_bin_bottom_0.5 more_cell_type_specific_mean_sad_bin_top_0.5 less_cell_type_specific_mean_sad_bin_bottom_0.3 less_cell_type_specific_mean_sad_bin_top_0.3 less_cell_type_specific_mean_sad_bin_bottom_0.1 less_cell_type_specific_mean_sad_bin_top_0.1 more_cell_type_specific_mean_sad_bin_bottom_0.3 more_cell_type_specific_mean_sad_bin_top_0.3 more_cell_type_specific_mean_sad_bin_bottom_0.1 more_cell_type_specific_mean_sad_bin_top_0.1
    
    # baseline without enhancer and promoter marks
    parallel -j 8 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --h2 ${enformer_dir}/ldsc/sumstats/${traits[$i]}.ldsc.imputed_v3.both_sexes.tsv.bgz --ref-ld-chr ${baseline_remove_enhancer_promoter_annot},${enformer_dir}/ldsc/annots/${tissues[$i]}_{}. --w-ld-chr ${weights} --frqfile-chr ${freqs} --overlap-annot --print-coefficients --out ${enformer_dir}/ldsc/out/${traits[$i]}_baseline_remove_enhancer_promoter_annots_${tissues[$i]}_{} ::: less_cell_type_specific less_cell_type_specific_mean_sad_bin_bottom_0.5 less_cell_type_specific_mean_sad_bin_top_0.5 more_cell_type_specific more_cell_type_specific_mean_sad_bin_bottom_0.5 more_cell_type_specific_mean_sad_bin_top_0.5 less_cell_type_specific_mean_sad_bin_bottom_0.3 less_cell_type_specific_mean_sad_bin_top_0.3 less_cell_type_specific_mean_sad_bin_bottom_0.1 less_cell_type_specific_mean_sad_bin_top_0.1 more_cell_type_specific_mean_sad_bin_bottom_0.3 more_cell_type_specific_mean_sad_bin_top_0.3 more_cell_type_specific_mean_sad_bin_bottom_0.1 more_cell_type_specific_mean_sad_bin_top_0.1
done

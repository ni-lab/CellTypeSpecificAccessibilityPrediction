#!/bin/bash

enformer_dir="/global/scratch/users/poojakathail/enformer"
freq_dir="/clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_frq"
ref_ld_list="/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/baselineLD.,${enformer_dir}/ldsc/annots/tissue_cell_type_specific_sad_by_tissue_tracks."

parallel -j 3 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --h2 ${enformer_dir}/ldsc/sumstats/{}.ldsc.imputed_v3.both_sexes.tsv.bgz --ref-ld-chr ${ref_ld_list} --w-ld-chr /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --frqfile-chr ${freq_dir}/1000G.EUR.QC. --overlap-annot --print-coefficients --out ${enformer_dir}/ldsc/out/{}_baseline_enformer_tissue_cell_type_specific_sad_by_tissue_tracks ::: 50_irnt 23104_irnt 6152_8 2443 6152_9 20116_0 78_irnt

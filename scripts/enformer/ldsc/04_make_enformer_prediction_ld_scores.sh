#!/bin/bash

enformer_dir="/global/scratch/users/poojakathail/enformer"

parallel -j 11 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --l2 --bfile /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink/1000G.EUR.QC.{} --ld-wind-cm 1 --annot ${enformer_dir}/ldsc/annots/tissue_cell_type_specific_sad_by_tissue_tracks_with_percentile_bins.{}.annot --thin-annot --out ${enformer_dir}/ldsc/annots/tissue_cell_type_specific_sad_by_tissue_tracks_with_percentile_bins..{} --print-snps /clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/list.txt ::: {1..22}
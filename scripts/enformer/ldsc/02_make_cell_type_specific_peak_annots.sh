#!/bin/bash

enformer_dir="/global/scratch/users/poojakathail/enformer"
g_dir="/clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink"
snp_list="/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots/list.txt"

for tissue in CNS Cardiovascular Blood_Immune Musculoskeletal_Connective Digestive Kidney Liver Pancreas;
do
	for bin in 0 1;
	do
		# liftover bed file
		/clusterfs/nilah/pooja/software/liftOver_linux ${enformer_dir}/ldsc/bed_files/${tissue}_bin_${bin}_cell_type_peak_sequences.hg38.bed /clusterfs/nilah/pooja/software/hg38ToHg19.over.chain ${enformer_dir}/ldsc/bed_files/${tissue}_bin_${bin}_cell_type_peak_sequences.hg19.bed ${enformer_dir}/ldsc/bed_files/${tissue}_bin_${bin}_cell_type_peak_sequences.hg38.bed.unmap

		# make annots
	    parallel -j 11 python /clusterfs/nilah/pooja/software/ldsc/make_annot.py \
	        --bed-file ${enformer_dir}/ldsc/bed_files/${tissue}_bin_${bin}_cell_type_peak_sequences.hg19.bed \
	        --bimfile ${g_dir}/1000G.EUR.QC.{}.bim \
	        --annot-file ${enformer_dir}/ldsc/annots/${tissue}_bin_${bin}_cell_type_peak_sequences.{}.annot.gz ::: {1..22}

	    # compute LD scores
	    parallel -j 11 python /clusterfs/nilah/pooja/software/ldsc/ldsc.py --l2 \
	    	--bfile ${g_dir}/1000G.EUR.QC.{} --ld-wind-cm 1 --thin-annot \
	    	--annot ${enformer_dir}/ldsc/annots/${tissue}_bin_${bin}_cell_type_peak_sequences.{}.annot.gz \ 
	    	--out ${enformer_dir}/ldsc/annots/${tissue}_bin_${bin}_cell_type_peak_sequences.{} --print-snps ${snp_list} ::: {1..22}
	done
done


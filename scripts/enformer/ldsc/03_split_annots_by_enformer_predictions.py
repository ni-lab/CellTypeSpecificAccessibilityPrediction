import numpy as np
import pandas as pd
import h5py

def main():
	enformer_dir = "/global/scratch/users/poojakathail/enformer"
	baseline_dir = "/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots"
	frq_dir = "/clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_frq"

	# mapping between traits and Enformer cell types
	enformer_targets = pd.read_csv(f"{enformer_dir}/targets_human_with_tissue_category.txt", sep="\t", index_col=0)
	dnase_atac_targets = enformer_targets[enformer_targets["description"].str.contains("DNASE") | 
                                      enformer_targets["description"].str.contains("ATAC")].index.values
	tissues = ['CNS',  'Cardiovascular', 'Blood/Immune', 'Musculoskeletal/Connective', 
	           'Digestive', 'Kidney', 'Liver', 'Pancreas']

	# read in enformer predictions for 1000 Genomes SNPs
	for chr in range(1, 23):
	    h5 = h5py.File(f"{enformer_dir}/variant-scores/1000-genomes/enformer/1000G.MAF_threshold=0.005.{chr}.h5", "r")
	    sad = h5["SAD"][:,dnase_atac_targets]
	    max_abs_sad = np.abs(sad).max(axis=1)
	    mean_abs_sad = np.abs(sad).mean(axis=1)
	    sad_df = pd.DataFrame({"chr": h5["chr"][:].astype(str),
	                           "pos": h5["pos"][:],
	                           "snp": h5["snp"][:].astype(str),
	                           "max_abs_sad": max_abs_sad,
	                           "mean_abs_sad": mean_abs_sad})
	    for tissue in tissues:
	        tissue_targets = enformer_targets[enformer_targets["Tissue category"] == tissue].index.values
	        sad_df[f"{tissue}_max_abs_sad"] = np.abs(sad[:, tissue_targets]).max(axis=1)
	        sad_df[f"{tissue}_mean_abs_sad"] = np.abs(sad[:, tissue_targets]).mean(axis=1)
	    sad_df.to_csv(f"{enformer_dir}/1000G.MAF_threshold=0.005.{chr}.sad_per_tissue.tsv", 
	                  sep="\t", header=True, index=False)

	# read in LDSC annotations of more and less cell type specific peak sequences for each tissue
	annots = {}
	for chr in range(1, 23):
	    baseline = pd.read_csv(f"{baseline_dir}/baselineLD.{chr}.annot.gz", sep="\t", compression="gzip")
	    sad = pd.read_csv(f"{enformer_dir}/1000G.MAF_threshold=0.005.{chr}.sad_per_tissue.tsv", sep="\t")
	    merged_df = baseline.merge(sad, left_on="SNP", right_on="snp", how="left")
	    
	    peak_annots = []
	    for tissue in tissues:
	        for bin in range(0, 2):
	        	if bin == 0:
	        		label = f"{tissue}_more_cell_type_specific"
	        	else:
	        		label = f"{tissue}_less_cell_type_specific"

	            bin_df = pd.read_csv(f"{enformer_dir}/ldsc/annots/{tissue}_bin_{bin}_cell_type_peak_sequences.{chr}.annot.gz",
	            					 compression="gzip", skiprows=1, names=[label])
	            peak_annots.append(bin_df)                           
	    peak_annots = pd.concat(peak_annots, axis=1)
	    merged_df = pd.concat([merged_df, peak_annots], axis=1)
	    annots[chr] = merged_df
	annots = pd.concat(annots.values())

	# split more and less cell type specific bins by Enformer predictions
	annot_cols = np.concatenate([[f"{tissue}_more_cell_type_specific", 
                              f"{tissue}_more_cell_type_specific_max_sad_bin_0",  
                              f"{tissue}_more_cell_type_specific_max_sad_bin_1",
                              f"{tissue}_more_cell_type_specific_mean_sad_bin_0",  
                              f"{tissue}_more_cell_type_specific_mean_sad_bin_1",
                              f"{tissue}_less_cell_type_specific", 
                              f"{tissue}_less_cell_type_specific_max_sad_bin_0", 
                              f"{tissue}_less_cell_type_specific_max_sad_bin_1",
                              f"{tissue}_less_cell_type_specific_mean_sad_bin_0", 
                              f"{tissue}_less_cell_type_specific_mean_sad_bin_1"] 
                             for tissue in tissues])

	for tissue in tissues:
		tissue_formatted = tissue.replace('_', '/')
	    for col in [f"{tissue}_more_cell_type_specific", f"{tissue}_less_cell_type_specific"]:
	        for sad_stat in ["max", "mean"]:
	        	# compute threshold of Enformer predictions to create two equally sized bins
	            bins = np.quantile(annots[annots[col] == 1][f"{tissue_formatted}_{sad_stat}_abs_sad"].values.flatten(),
	                               np.arange(0, 1.5, 0.5))
	            for i in range(len(bins) - 1):
	                annots[f"{col}_{sad_stat}_sad_bin_{i}"] = 0
	                annots[f"{col}_{sad_stat}_sad_bin_{i}"].loc[(annots[col] == 1) &
	                                                 			(annots[f"{tissue_formatted}_{sad_stat}_abs_sad"] >= bins[i]) &
	                                                 			(annots[f"{tissue_formatted}_{sad_stat}_abs_sad"] < bins[i+1])] = 1
	# write out annotations and frequnecy files per chromosome for LDSC
	for chr in range(1, 23):
	    chr_annots = annots[annots["CHR"] == chr]
	    chr_annots.index = np.arange(len(chr_annots))
	    chr_annots[annot_cols].to_csv(f"{enformer_dir}/ldsc/annots/tissue_cell_type_specific_sad_by_tissue_tracks.{chr}.annot",
	                                  sep="\t", header=True, index=False)

	    ## compute M and M_5_50 files
	    frq_df = pd.read_csv(f"{frq_dir}/1000G.EUR.QC.{chr}.frq",
	                         sep='\s+', header=0, index_col=False)

	    assert np.all(frq_df["SNP"].values == chr_annots["SNP"].values)

	    with open(f"{enformer_dir}/ldsc/annots/tissue_cell_type_specific_sad_by_tissue_tracks.{chr}.l2.M", "w") as f:
	        f.write("\t".join(chr_annots[annot_cols].sum().astype(str).values))


	    with open(f"{enformer_dir}/ldsc/annots/tissue_cell_type_specific_sad_by_tissue_tracks.{chr}.l2.M_5_50", "w") as f:
	        f.write("\t".join(chr_annots[annot_cols][frq_df["MAF"] > 0.05].sum().astype(str).values))


################################################################################
if __name__ == '__main__':
	main()


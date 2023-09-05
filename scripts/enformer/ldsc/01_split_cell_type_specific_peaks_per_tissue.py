import numpy as np
import pandas as pd

def main():
	enformer_dir = "/global/scratch/users/poojakathail/enformer"

	# mapping between traits and Enformer cell types
	enformer_targets = pd.read_csv(f"{enformer_dir}/targets_human_with_tissue_category.txt", sep="\t", index_col=0)

	# read in sequence_bins 
	seq_targets = pd.read_csv(f"{enformer_dir}/all_sequence_bins.csv", sep="\t", 
	                          names=["Chromosome", "Start", "End", "Name"])

	## read in processed target data
	with open(f"{enformer_dir}/train_targets_peaks.npz", "rb") as f:
	    train_targets_peaks = np.load(f)
	with open(f"{enformer_dir}/val_test_targets_peaks.npz", "rb") as f:
	    val_test_targets_peaks = np.load(f) 
	all_targets_peaks = np.concatenate([train_targets_peaks, val_test_targets_peaks])
	del train_targets_peaks, val_test_targets_peaks
	peak_in_num_cell_types = all_targets_peaks.sum(axis=1)

	tissues = ['CNS',  'Cardiovascular', 'Blood/Immune', 'Musculoskeletal/Connective', 
	           'Digestive', 'Kidney', 'Liver', 'Pancreas']
	proportion_tissue_peaks_threshold = 0.3

	# write bed files for more and less cell type specific peak sequences for each tissue
	for i, tissue in enumerate(tissues):
	    tissue_formatted = tissue.replace("/", "_")
	    tissue_inds = enformer_targets[enformer_targets["Tissue category"] == tissue].index.values
	    tissue_peaks = all_targets_peaks[:, tissue_inds].sum(axis=1)/len(tissue_inds)
	    tissue_peak_inds = np.where((tissue_peaks > proportion_tissue_peaks_threshold))[0]
	    
	    bins = np.percentile(peak_in_num_cell_types[tissue_peak_inds], np.arange(0, 150, 50))
	    for i in range(len(bins) - 1):
	        if i == len(bins) - 2:
	            bin_inds = np.where((peak_in_num_cell_types >= bins[i]) & 
	                                (peak_in_num_cell_types <= bins[i+1]))[0]
	        else:
	            bin_inds = np.where((peak_in_num_cell_types >= bins[i]) & 
	                                (peak_in_num_cell_types < bins[i+1]))[0]
	        bin_inds = np.intersect1d(bin_inds, tissue_peak_inds)
	        seq_targets.iloc[bin_inds].to_csv(f"{enformer_dir}/ldsc/bed_files/{tissue_formatted}_bin_{i}_cell_type_peak_sequences.hg38.bed", sep="\t", header=False, index=False)


################################################################################
if __name__ == '__main__':
	main()


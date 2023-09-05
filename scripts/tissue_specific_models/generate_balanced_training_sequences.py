from optparse import OptionParser
import os
import numpy as np
import pandas as pd

'''
generate_balanced_training_sequences.py 

Generate a balanced training set of peaks and non-peaks.
'''
def main():
	usage = 'usage: %prog [options]'
	parser = OptionParser(usage)
	parser.add_option('-p', dest='peaks',
	    help='Sequences overlapping peaks.')
	parser.add_option('-n', dest='non_peaks',
	    help='Sequences not overlapping peaks.')
	parser.add_option('-o', dest='out_dir',
	    help='Output directory.')
	parser.add_option('--gc', dest='gc_matched',
	    default=False, action='store_true', help='Match GC content of non-peaks to peaks.')
	(options, args) = parser.parse_args()

	NUM_BINS=20

	all_peaks = pd.read_csv(options.peaks, sep="\t")
	non_peaks = pd.read_csv(options.non_peaks, sep="\t")

	# sample non-peak sequences per fold
	nonpeaks_sampled = []
	for fold in ["train", "valid", "test"]:
	    sample_size = all_peaks[all_peaks["4_usercol"] == fold].shape[0]
	    non_peak_fold = non_peaks[non_peaks["4_usercol"] == fold]

	    if options.gc_matched:
	    	# assign sequences to GC content bins
	    	peak_gc = all_peaks[all_peaks["4_usercol"] == fold]["6_pct_gc"].values
		    non_peak_gc = non_peaks[non_peaks["4_usercol"] == fold]["6_pct_gc"].values

		    bins = np.linspace(np.min(peak_gc), np.max(peak_gc), num=NUM_BINS)
		    peak_binned = np.digitize(peak_gc, bins=bins)
		    nonpeak_binned = np.digitize(non_peak_gc, bins=bins)
		    peak_bin_is, peak_counts = np.unique(peak_binned, return_counts=True)
		    nonpeak_bin_is, nonpeak_counts = np.unique(nonpeak_binned, return_counts=True)
		    match_pos = np.isin(nonpeak_bin_is, peak_bin_is)
		    nonpeak_counts = nonpeak_counts[match_pos]

		    # Sample from bins
		    sample_size = len(peak_gc)
		    rng = np.random.default_rng(seed=0)
		    gc_matched_nonpeaks_fold = []
		    for bin_i in peak_bin_is:
		        peak_set_i = peak_gc[peak_binned == bin_i]
		        nonpeak_set_i = np.where(nonpeak_binned == bin_i)[0]
		        if nonpeak_set_i.shape[0] > 0:
		            if nonpeak_set_i.shape[0] >= peak_set_i.shape[0]:
		                nonpeak_set_i = rng.choice(nonpeak_set_i, size=peak_set_i.shape[0],
		                                           replace=False, axis=0)
		            else:
		                nonpeak_set_i = rng.choice(nonpeak_set_i, size=peak_set_i.shape[0],
		                                           replace=True, axis=0)
		            gc_matched_nonpeaks_fold.append(nonpeak_set_i)
		        else:
		            print(fold, bin_i)
		    gc_matched_nonpeaks_fold = np.concatenate(gc_matched_nonpeaks_fold, axis=0)
		    nonpeaks_sampled.append(non_peak_fold.iloc[gc_matched_nonpeaks_fold])

	    else:
	    	nonpeaks_sampled.append(non_peak_fold.sample(sample_size, replace=False))

	nonpeaks_sampled = pd.concat(nonpeaks_sampled)
	 
	# shuffle sequences within each fold
	training_seqs = pd.concat([all_peaks, nonpeaks_sampled])
	training_seqs_shuffled = []
	for fold in ["train", "valid", "test"]:
	    training_seqs_shuffled.append(training_seqs[training_seqs["4_usercol"] == fold].sample(frac=1, replace=False))
	training_seqs_shuffled = pd.concat(training_seqs_shuffled)

	# write out sequences
	if not os.path.isdir(options.out_dir):
    	os.mkdir(options.out_dir)

	training_seqs_shuffled[["#1_usercol", "2_usercol", "3_usercol", "4_usercol"]].to_csv(
	    f"{options.out_dir}/sequences.bed",
	    sep="\t", header=False, index=False)


################################################################################
if __name__ == '__main__':
	main()


import os
import numpy as np
import pandas as pd
import h5py

"""
After running 04_make_enformer_prediction_ld_scores.sh, run this script to split the LD scores for each annotation 
into a separate file so that they can be run independently.
"""

def main():
    annot_dir = "/global/scratch/users/poojakathail/enformer/ldsc/annots"
    annot_prefix = "tissue_cell_type_specific_sad_by_tissue_tracks_with_percentile_bins"
    
    for chr in range(1, 23):
        print("Writing LD scores for chr", chr)
        
        chr_ld_scores = pd.read_csv(f"{annot_dir}/{annot_prefix}.{chr}.l2.ldscore.gz", 
                                    sep="\t", compression="gzip", index_col=[0,1,2])
        m = pd.read_csv(f"{annot_dir}/{annot_prefix}.{chr}.l2.M", 
                                    sep="\t", header=None)
        m_5_50 = pd.read_csv(f"{annot_dir}/{annot_prefix}.{chr}.l2.M_5_50", 
                                    sep="\t", header=None)
        
        for i, col in enumerate(chr_ld_scores.columns):
            col_file_name = col[:-2].replace("/", "_")
            chr_ld_scores[col].to_csv(f"{annot_dir}/{col_file_name}.{chr}.l2.ldscore.gz",
                                      sep="\t", compression="gzip", header=True, index=True)
            m[i].to_csv(f"{annot_dir}/{col_file_name}.{chr}.l2.M",
                          sep="\t", header=False, index=False)
            m_5_50[i].to_csv(f"{annot_dir}/{col_file_name}.{chr}.l2.M_5_50",
                          sep="\t", header=False, index=False)
            
################################################################################
if __name__ == '__main__':
    main()


import os
import numpy as np
import pandas as pd
import h5py
from tqdm import tqdm

from optparse import OptionParser

def main():
    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)
    parser.add_option('--tissue', dest='tissue',
        help='Tissue.')
    parser.add_option('--dataset', dest='dataset',
        help='Dataset.')
    parser.add_option('--tasks', dest='tasks',
        help='Tasks.')
    parser.add_option('--targets', dest='targets',
        help='targets.')
    (options, args) = parser.parse_args()
    
    enformer_dir = "/global/scratch/users/poojakathail/enformer"
    baseline_dir = "/clusterfs/nilah/pooja/immune_atlas/ldsc/baseline_annots"
    frq_dir = "/clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_Phase3_frq"
    tissue = options.tissue
    dataset = options.dataset
    tasks = options.tasks
    targets = options.targets.split(",")
    targets_df = pd.read_csv(f"/clusterfs/nilah/pooja/{dataset}/targets.txt", sep="\t", index_col=0)
    targets_i = targets_df[targets_df["identifier"].isin(targets)].index.values
    
    print("Processing tissue specific predictions for 1000 genomes SNPs..")
    annots = {}
    for chr in tqdm(range(1, 23)):
        baseline = pd.read_csv(f"{baseline_dir}/baselineLD.{chr}.annot.gz", sep="\t", compression="gzip")
        
        if tasks in ["multitask", "multitask_8x_params_same_layers"]:
            h5 = h5py.File(
                f"/clusterfs/nilah/pooja/{dataset}/ldsc/sad/train__{tasks}__all_sequences__1/chr{chr}/sad.h5", 
                "r")
            sad = h5["SAD"][:,targets_i]
            assert ~np.any(np.isnan(sad))
            assert ~np.any(np.isinf(sad))
            # if np.any(np.isnan(sad)):
            #     print(chr)
            # sad = np.nan_to_num(sad)
            
            max_abs_sad = np.abs(sad).max(axis=1)
            mean_abs_sad = np.abs(sad).mean(axis=1)
            sad_df = pd.DataFrame({"chr": h5["chr"][:].astype(str),
                                   "pos": h5["pos"][:],
                                   "snp": h5["snp"][:].astype(str),
                                   "max_abs_sad": max_abs_sad,
                                   "mean_abs_sad": mean_abs_sad})
        else:
            sad_df = pd.DataFrame([])
            for i, cell_type in enumerate(targets):
                h5 = h5py.File(
                      f"/clusterfs/nilah/pooja/{dataset}/ldsc/sad/train__{tasks}__{cell_type}__all_sequences__1/chr{chr}/sad.h5", 
                      "r")
                sad = h5["SAD"][:,0]
                assert ~np.any(np.isnan(sad))
                assert ~np.any(np.isinf(sad))
                # if np.any(np.isnan(sad)):
                #     print(cell_type, chr)
                # sad = np.nan_to_num(sad)
                
                if "snp" not in sad_df.columns:
                    sad_df["chr"] = h5["chr"][:].astype(str)
                    sad_df["pos"] = h5["pos"][:]
                    sad_df["snp"] = h5["snp"][:].astype(str)

                assert np.all(sad_df["snp"].values == h5["snp"][:].astype(str))
                sad_df[f"{cell_type}_sad"] = sad
            sad_df["max_abs_sad"] = sad_df[[f"{cell_type}_sad" for cell_type in targets]].abs().max(axis=1)
            sad_df["mean_abs_sad"] = sad_df[[f"{cell_type}_sad" for cell_type in targets]].abs().mean(axis=1)
        
        merged_df = baseline.merge(sad_df, left_on="SNP", right_on="snp", how="left")
        merged_df = merged_df.fillna(0) # if model is trained in hg38, some SNPs could not be lifted over so predictions will be NaN
        
        # read in LDSC annotations of more and less cell type specific peak sequences for each tissue
        peak_annots = []
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

    print("Split more and less cell type specific bins by tissue specific predictions..")
    annot_cols = []
    for col in [f"{tissue}_more_cell_type_specific", f"{tissue}_less_cell_type_specific"]:
        for sad_stat in ["max", "mean"]:
            for percentile in [0.5]:
                percentile_bins = np.quantile(annots[annots[col] == 1][f"{sad_stat}_abs_sad"].values.flatten(),
                                                    [percentile, 1-percentile])
                print(annots[col].sum())
                print(percentile_bins)
                
                annots[f"{col}_{dataset}_{tasks}_{sad_stat}_sad_bin_bottom_{percentile}"] = 0
                annots[f"{col}_{dataset}_{tasks}_{sad_stat}_sad_bin_bottom_{percentile}"].loc[(annots[col] == 1) &
                                                            (annots[f"{sad_stat}_abs_sad"] <= percentile_bins[0])] = 1

                annots[f"{col}_{dataset}_{tasks}_{sad_stat}_sad_bin_top_{percentile}"] = 0
                annots[f"{col}_{dataset}_{tasks}_{sad_stat}_sad_bin_top_{percentile}"].loc[(annots[col] == 1) &
                                                            (annots[f"{sad_stat}_abs_sad"] >= percentile_bins[1])] = 1
                annot_cols.extend([f"{col}_{dataset}_{tasks}_{sad_stat}_sad_bin_bottom_{percentile}", 
                                   f"{col}_{dataset}_{tasks}_{sad_stat}_sad_bin_top_{percentile}"])

    print(f"Writing annotations and frequency files per chromosome for LDSC..")
    for chr in tqdm(range(1, 23)):
        chr_annots = annots[annots["CHR"] == chr]
        chr_annots.index = np.arange(len(chr_annots))

        ## compute M and M_5_50 files
        frq_df = pd.read_csv(f"{frq_dir}/1000G.EUR.QC.{chr}.frq",
                             sep='\s+', header=0, index_col=False)
        assert np.all(frq_df["SNP"].values == chr_annots["SNP"].values)

        for annot_col in annot_cols:
            chr_annots[[annot_col]].to_csv(f"/clusterfs/nilah/pooja/{dataset}/ldsc/annots/{annot_col}.{chr}.annot",
                                        sep="\t", header=True, index=False)

            with open(f"/clusterfs/nilah/pooja/{dataset}/ldsc/annots/{annot_col}.{chr}.l2.M", "w") as f:
                f.write("\t".join(chr_annots[annot_cols].sum().astype(str).values))


            with open(f"/clusterfs/nilah/pooja/{dataset}/ldsc/annots/{annot_col}.{chr}.l2.M_5_50", "w") as f:
                f.write("\t".join(chr_annots[annot_cols][frq_df["MAF"] > 0.05].sum().astype(str).values))


################################################################################
if __name__ == '__main__':
    main()


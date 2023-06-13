import h5py
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn.ensemble import RandomForestClassifier
from tqdm import tqdm

enformer_dir = "/global/scratch/users/poojakathail/enformer"
bins = [(1, 10), (10, 50), (50, 100),
        (100, 300), (300, 685)]

bin_edges = [1, 10, 50, 100, 300, 685]

def read_enformer_predictions(tissues):
    pos_preds = {}
    neg_preds = {}
    for tissue in tissues:
        try:
            pos_preds[tissue] = h5py.File(f"{enformer_dir}/gtex_fine/preds/{tissue}_pos_REF_ALT/sad.h5", "r")
            neg_preds[tissue] = h5py.File(f"{enformer_dir}/gtex_fine/preds/{tissue}_neg_REF_ALT/sad.h5", "r")
        except Exception as e:
            continue
    return {"pos": pos_preds,
            "neg": neg_preds}
        
def read_sei_predictions(tissues):  
    pos_preds = {}
    neg_preds = {}
    for tissue in tissues:
        try:
            pos_preds[tissue] = h5py.File(f"{enformer_dir}/gtex_fine/preds/sei_{tissue}_pos/chromatin-profiles-hdf5/{tissue}_pos_diffs.h5", "r")
            neg_preds[tissue] = h5py.File(f"{enformer_dir}/gtex_fine/preds/sei_{tissue}_neg/chromatin-profiles-hdf5/{tissue}_neg_diffs.h5", "r")
            if pos_preds[tissue]["data"][:,:].sum() == 0 or neg_preds[tissue]["data"][:,:].sum() == 0:
                del pos_preds[tissue]
                del neg_preds[tissue]
        except Exception as e:
            continue
    return {"pos": pos_preds,
            "neg": neg_preds}


def stratified_classification_performance_all_tissues(tissues, enformer_preds, sei_preds,
                                                      dnase_atac_targets, peak_in_num_cell_types):

    classification_performance = []
    for tissue in tqdm(tissues):
        if tissue in enformer_preds["pos"] and tissue in enformer_preds["neg"] and tissue in sei_preds["pos"] and tissue in sei_preds["neg"]:
            classification_performance.append(stratified_classification_performance_for_tissue(tissue, enformer_preds, 
                                                                                               sei_preds,
                                                                                               dnase_atac_targets, 
                                                                                               peak_in_num_cell_types))
    classification_performance = pd.concat(classification_performance)
    return classification_performance

def stratified_classification_performance_for_tissue(tissue, enformer_preds, sei_preds, 
                                                     dnase_atac_targets, peak_in_num_cell_types):
    # read in data mapping vcf to basenji bins
    pos_set = pd.read_csv(f"{enformer_dir}/gtex_fine/vcf/{tissue}_pos_basenji_bins_intersect.tsv",
                        sep="\t", index_col=False, 
                        names=["Chromosome", "Start", "snp", "a1", "a2", ".", "..", "na", 
                              "Bin Chromosome", "Bin Start", "Bin End", "Bin Name"])
    pos_set = pos_set.drop_duplicates(subset=["Chromosome", "Start", "snp", "a1", "a2"])
    pos_set = pos_set.reset_index()
    
    # join to get beta sign and TSS distance
    susie_results = pd.read_csv(f"{enformer_dir}/gtex_fine/susie/{tissue}.tsv", sep="\t", index_col=0)
    gtex_genes = pd.read_csv("/clusterfs/nilah/pooja/genomes/gtex_genes.bed",
                            sep="\t", index_col=0)[["name", "chrom", "txStart"]]
    pos_set["index"] = pos_set.index.values
    pos_set = pos_set.merge(susie_results, left_on="snp", right_on="variant", how="left")
    pos_set = pos_set.merge(gtex_genes, left_on="gene", right_on="name")
    pos_set["TSS distance"] = (pos_set["txStart"] - pos_set["Start"]).abs()

    # filter to only eqtls with pip > 0.9
    pos_set = pos_set[pos_set["pip"] > 0.9]

    # retain closest gene per variant
    pos_set = pos_set.loc[pos_set.groupby('variant')["TSS distance"].idxmin()]
    pos_set.index = pos_set["index"]
                            
    neg_set = pd.read_csv(f"{enformer_dir}/gtex_fine/vcf/{tissue}_neg_basenji_bins_intersect.tsv",
                        sep="\t", index_col=False, 
                        names=["Chromosome", "Start", "snp", "a1", "a2", ".", "..", "na",
                              "Bin Chromosome", "Bin Start", "Bin End", "Bin Name"])
    neg_set = neg_set.drop_duplicates(subset=["Chromosome", "Start", "snp", "a1", "a2"])
    neg_set = neg_set.reset_index()

    # subset to variants in a peak
    pos_set_matches = pos_set[pos_set["Bin Name"] != "."].index.values
    pos_set_matches = pos_set_matches[np.where(peak_in_num_cell_types[pos_set.loc[pos_set_matches]["Bin Name"].values.astype(int)] > 0)[0]]

    # split variants into bins based on cell type specificity of peak
    bin_assignments = np.digitize(peak_in_num_cell_types[pos_set.loc[pos_set_matches]["Bin Name"].values.astype(int)], bin_edges)
    
    classification_performance = pd.DataFrame([])
    if len(pos_set_matches) > 500:
        # compute performance in each bin
        for bin_i in range(1, len(bins)):
            pos_set_inds = np.sort(pos_set_matches[np.where(bin_assignments == bin_i)[0]])
            neg_set_inds = np.sort(np.random.choice(len(neg_set), size=len(pos_set_inds), replace=False))

            ## enformer
            X = np.vstack([enformer_preds["pos"][tissue]["SAD"][pos_set_inds,:], 
                           enformer_preds["neg"][tissue]["SAD"][neg_set_inds,:]])
            X = X[:,dnase_atac_targets["enformer"]]
            y = np.concatenate([np.ones(len(pos_set_inds)),
                                np.zeros(len(neg_set_inds))])

            rf = RandomForestClassifier(max_features="log2")
            enformer_cv_scores = cross_validate(rf, X, y, scoring=["roc_auc", "average_precision"], cv=8, n_jobs=8)
            classification_performance = classification_performance.append({"Tissue": tissue,
                                                                            "Bin": bins[bin_i-1],
                                                                            "Model": "Enformer",
                                                                            "AUC": np.mean(enformer_cv_scores["test_roc_auc"]),
                                                                            "AUPRC": np.mean(enformer_cv_scores["test_average_precision"]),
                                                                           "Positive set size": len(pos_set_inds),
                                                                            "TSS distance": list(pos_set.loc[pos_set_inds, "TSS distance"].values),
                                                                           "Beta": list(pos_set.loc[pos_set_inds, "beta_marginal"].values)},
                                                                            ignore_index=True)

            ## sei
            X = np.vstack([sei_preds["pos"][tissue]["data"][pos_set_inds,:], 
                           sei_preds["neg"][tissue]["data"][neg_set_inds,:]])
            X = X[:,dnase_atac_targets["sei"]]
            y = np.concatenate([np.ones(len(pos_set_inds)),
                                np.zeros(len(neg_set_inds))])

            rf = RandomForestClassifier(max_features="log2")
            sei_cv_scores = cross_validate(rf, X, y, scoring=["roc_auc", "average_precision"], cv=8, n_jobs=8)
            classification_performance = classification_performance.append({"Tissue": tissue,
                                                                            "Bin": bins[bin_i-1],
                                                                            "Model": "Sei",
                                                                            "AUC": np.mean(sei_cv_scores["test_roc_auc"]),
                                                                            "AUPRC": np.mean(sei_cv_scores["test_average_precision"]),
                                                                           "Positive set size": len(pos_set_inds),
                                                                            "TSS distance": list(pos_set.loc[pos_set_inds, "TSS distance"].values),
                                                                           "Beta": list(pos_set.loc[pos_set_inds, "beta_marginal"].values)},
                                                                            ignore_index=True)
    return classification_performance        

        
def stratified_direction_performance_all_tissues(tissues, enformer_preds, sei_preds,
                                                      dnase_atac_targets, peak_in_num_cell_types):

    direction_performance = []
    for tissue in tqdm(tissues):
        if tissue in enformer_preds["pos"] and tissue in sei_preds["pos"]:
            direction_performance.append(stratified_direction_performance_for_tissue(tissue, enformer_preds, 
                                                                                               sei_preds,
                                                                                               dnase_atac_targets, 
                                                                                               peak_in_num_cell_types))
    direction_performance = pd.concat(direction_performance)
    return direction_performance



        
def stratified_direction_performance_for_tissue(tissue, enformer_preds, sei_preds, 
                                                     dnase_atac_targets, peak_in_num_cell_types):
    # read in data mapping vcf to basenji bins
    pos_set = pd.read_csv(f"{enformer_dir}/gtex_fine/vcf/{tissue}_pos_basenji_bins_intersect.tsv",
                        sep="\t", index_col=False, 
                        names=["Chromosome", "Start", "snp", "a1", "a2", ".", "..", "na", 
                              "Bin Chromosome", "Bin Start", "Bin End", "Bin Name"])
    pos_set = pos_set.drop_duplicates(subset=["Chromosome", "Start", "snp", "a1", "a2"])
    pos_set = pos_set.reset_index()
    
    # join to get beta sign and TSS distance
    susie_results = pd.read_csv(f"{enformer_dir}/gtex_fine/susie/{tissue}.tsv", sep="\t", index_col=0)
    gtex_genes = pd.read_csv("/clusterfs/nilah/pooja/genomes/gtex_genes.bed",
                            sep="\t", index_col=0)[["name", "chrom", "txStart"]]
    pos_set["index"] = pos_set.index.values
    pos_set = pos_set.merge(susie_results, left_on="snp", right_on="variant", how="left")
    pos_set = pos_set.merge(gtex_genes, left_on="gene", right_on="name")
    pos_set["eQTL direction"] = pos_set["beta_marginal"] > 0
    pos_set["TSS distance"] = (pos_set["txStart"] - pos_set["Start"]).abs()

    # filter to only variant that have the same sign across multiple cis genes
    variant_mask = pos_set.groupby('variant')["eQTL direction"].unique().apply(lambda x: len(x))
    variant_mask = variant_mask[variant_mask == 1].index.values
    pos_set = pos_set[pos_set["variant"].isin(variant_mask)]

    # filter to only eqtls with pip > 0.9
    pos_set = pos_set[pos_set["pip"] > 0.9]

    # retain closest gene per variant
    pos_set = pos_set.loc[pos_set.groupby('variant')["TSS distance"].idxmin()]
    pos_set.index = pos_set["index"]


    # subset to variants in a peak
    pos_set_matches = pos_set[pos_set["Bin Name"] != "."].index.values
    pos_set_matches = pos_set_matches[np.where(peak_in_num_cell_types[pos_set.loc[pos_set_matches]["Bin Name"].values.astype(int)] > 0)[0]]

    # split variants into bins based on cell type specificity of peak
    bin_assignments = np.digitize(peak_in_num_cell_types[pos_set.loc[pos_set_matches]["Bin Name"].values.astype(int)], bin_edges)
    
    direction_performance = pd.DataFrame([])
    if len(pos_set_matches) > 1000:
        # compute performance in each bin
        for bin_i in range(1, len(bins)):
            pos_set_inds = np.sort(pos_set_matches[np.where(bin_assignments == bin_i)[0]])

            ## enformer
            X = enformer_preds["pos"][tissue]["SAD"][pos_set_inds,:]
            X = X[:,dnase_atac_targets["enformer"]]
            y = pos_set.loc[pos_set_inds, "eQTL direction"]

            rf = RandomForestClassifier(max_features="log2")
            enformer_cv_scores = cross_validate(rf, X, y, scoring=["roc_auc", "average_precision"], cv=8, n_jobs=8)
            direction_performance = direction_performance.append({"Tissue": tissue,
                                                                            "Bin": bins[bin_i-1],
                                                                            "Model": "Enformer",
                                                                            "AUC": np.mean(enformer_cv_scores["test_roc_auc"]),
                                                                            "AUPRC": np.mean(enformer_cv_scores["test_average_precision"]),
                                                                 "Test set size": len(pos_set_inds)},                                                                        ignore_index=True)

            ## sei
            X = sei_preds["pos"][tissue]["data"][pos_set_inds,:]
            X = X[:,dnase_atac_targets["sei"]]
            y = pos_set.loc[pos_set_inds, "eQTL direction"]

            rf = RandomForestClassifier(max_features="log2")
            sei_cv_scores = cross_validate(rf, X, y, scoring=["roc_auc", "average_precision"], cv=8, n_jobs=8)
            direction_performance = direction_performance.append({"Tissue": tissue,
                                                                            "Bin": bins[bin_i-1],
                                                                            "Model": "Sei",
                                                                            "AUC": np.mean(sei_cv_scores["test_roc_auc"]),
                                                                            "AUPRC": np.mean(sei_cv_scores["test_average_precision"]),
                                                                 "Test set size": len(pos_set_inds)},
                                                                            ignore_index=True)
    return direction_performance 

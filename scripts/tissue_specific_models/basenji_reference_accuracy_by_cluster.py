import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import h5py

from optparse import OptionParser


def main():
    usage = 'usage: %prog [options] <output_file> <targets_dir>'
    parser = OptionParser(usage)
    parser.add_option('-t', dest='targets_file',
        help='Targets file.')
    
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide output_file and targets directory.')
    else:
        output_file = args[0]
        model_dir = output_file.split("/predict_beds")[0]
        targets_dir = args[1]

    targets = pd.read_csv(options.targets_file, sep="\t", header=0, index_col=0)
    cell_types = targets["identifier"].values
    clusters = [d for d in os.listdir(f"{model_dir}/predict_beds") if os.path.isdir(f"{model_dir}/predict_beds/{d}")]

    log_performance_by_cluster = pd.DataFrame([], index=clusters, columns=cell_types)
    for cluster in clusters:
        preds = h5py.File(f"{model_dir}/predict_beds/{cluster}/predict.h5", "r")
        preds_df = pd.DataFrame(np.nan_to_num(np.squeeze(preds["preds"][:,:,:])), columns=[f"{ct}_pred" for ct in cell_types])  
        preds_df["chrom"] = preds["chrom"][:].astype(str)
        preds_df["start"] = preds["start"][:]
        preds_df["end"] = preds["end"][:]

        predict_regions = pd.read_csv(f"{targets_dir}/{cluster}/predict_regions.bed", sep="\t", names=["chrom", "start", "end", "name"])
        preds_df = preds_df.merge(predict_regions, on=["chrom", "start", "end"], how="inner")

        for cell_type in cell_types:
            cell_type_targets = pd.read_csv(f"{targets_dir}/{cluster}/{cell_type}_target_signal.out", sep="\t", 
                                            index_col=0, names=["size", "covered", "sum", "mean0", "mean"])
            cell_type_targets = cell_type_targets.loc[preds_df["name"].values]["sum"].values
            preds_df[f"{cell_type}_target"] = cell_type_targets

            log_performance_by_cluster.loc[cluster, cell_type] = pearsonr(np.log2(preds_df[f"{cell_type}_target"]+1),
                                                                          np.log2(preds_df[f"{cell_type}_pred"]+1))[0]


    log_performance_by_cluster.to_csv(output_file, sep="\t", header=True, index=True)
        


if __name__ == '__main__':
  main()






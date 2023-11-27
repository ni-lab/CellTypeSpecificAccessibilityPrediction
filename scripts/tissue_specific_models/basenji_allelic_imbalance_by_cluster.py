import os
import numpy as np
import pandas as pd
import h5py

from scipy.stats import pearsonr
from sklearn.metrics import roc_curve
from compare_auc_delong_xu import delong_roc_variance

from optparse import OptionParser


def main():
    usage = 'usage: %prog [options] <output_file>'
    parser = OptionParser(usage)
    parser.add_option('-t', dest='targets_file',
        help='Targets file.')
    parser.add_option('-i', dest='imbalance_dir',
        help='Allelic imbalance directory.')
    
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide output_file.')
    else:
        output_file = args[0]
        model_dir = output_file.split("/allelic_imbalance")[0]

    targets = pd.read_csv(options.targets_file, sep="\t", header=0, index_col=0)
    cell_types = [d for d in os.listdir(f"{model_dir}/allelic_imbalance") if os.path.isdir(f"{model_dir}/allelic_imbalance/{d}")]
    vcf_columns = ['chr', 'pos', 'snp', 'Ref', 'Alt', '#Ref', '#Alt', 'P-value', 'imbalance', '#Total', 'Q-value']

    allelic_imbalance_auc_by_cluster = pd.DataFrame([], columns=cell_types)
    allelic_imbalance_cov_by_cluster = pd.DataFrame([], columns=cell_types)
    for cell_type in cell_types:
        ti = targets[targets["identifier"] == cell_type].index.values[0]
        clusters = [d for d in os.listdir(f"{model_dir}/allelic_imbalance/{cell_type}") if os.path.isdir(f"{model_dir}/allelic_imbalance/{cell_type}/{d}")]

        for cluster in clusters:
            vcf = pd.read_csv(f"{options.imbalance_dir}/{cell_type}/{cluster}.vcf", sep="\t", names=vcf_columns)
            preds = h5py.File(f"{model_dir}/allelic_imbalance/{cell_type}/{cluster}/sad.h5", "r")
            ref_preds = preds["REF"][:,:,ti].flatten().astype(np.float32)
            alt_preds = preds["ALT"][:,:,ti].flatten().astype(np.float32)

            imb_labels = (vcf["imbalance"] > 0.5).values
            imb_preds = ref_preds/(ref_preds+alt_preds)
            fpr, tpr, thresholds = roc_curve(imb_labels, imb_preds)
            roc_auc, auc_cov = delong_roc_variance(imb_labels, imb_preds)
            allelic_imbalance_auc_by_cluster.loc[cluster, cell_type] = roc_auc
            allelic_imbalance_cov_by_cluster.loc[cluster, cell_type] = auc_cov

    with pd.ExcelWriter(output_file) as writer:  
        allelic_imbalance_auc_by_cluster.to_excel(writer, sheet_name='roc_auc')
        allelic_imbalance_cov_by_cluster.to_excel(writer, sheet_name='auc_cov')


if __name__ == '__main__':
  main()






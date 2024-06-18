import os
import numpy as np
import pandas as pd
import h5py
from scipy.stats import pearsonr


def read_reference_accuracy_replicates(data_dir, cell_types, tasks, training_sets):
    ref_acc_by_cluster = {}
    for task in tasks:
        ref_acc_by_cluster[task] = {}
        for training_set in training_sets:
            ref_acc_by_cluster[task][training_set] = []
            for replicate in range(1, 4):
                if "single" in task:
                    rep_data = []
                    for cell_type in cell_types:
                        try:
                            rep_data.append(pd.read_csv(f"{data_dir}/train/replicate_models/train__{task}__{cell_type}__{training_set}__{replicate}/predict_beds/reference_accuracy_by_cluster.tsv", 
                                                                         sep="\t", index_col=0))
                        except:
                            continue
                    if len(rep_data) > 0:   
                        rep_data = pd.concat(rep_data, axis=1)
                    ref_acc_by_cluster[task][training_set].append(rep_data)
                else:
                    try:
                        ref_acc_by_cluster[task][training_set].append(pd.read_csv(f"{data_dir}/train/replicate_models/train__{task}__{training_set}__{replicate}/predict_beds/reference_accuracy_by_cluster.tsv", 
                                                                             sep="\t", index_col=0))
                    except:
                        ref_acc_by_cluster[task][training_set].append([])
    return ref_acc_by_cluster


def read_imbalance_accuracy_replicates(data_dir, cell_types, tasks, training_sets, file_name="allelic_imbalance_by_cluster.xlsx"):
    ref_acc_by_cluster = {}
    for task in tasks:
        ref_acc_by_cluster[task] = {}
        for training_set in training_sets:
            ref_acc_by_cluster[task][training_set] = []
            for replicate in range(1, 4):
                if "single" in task:
                    rep_data = []
                    for cell_type in cell_types:
                        try:
                            rep_data.append(pd.read_excel(f"{data_dir}/train/replicate_models/train__{task}__{cell_type}__{training_set}__{replicate}/allelic_imbalance/{file_name}", 
                                                                         sheet_name="roc_auc", index_col=0)[cell_type])
                        except:
                            continue
                    if len(rep_data) > 0:   
                        rep_data = pd.concat(rep_data, axis=1)
                    ref_acc_by_cluster[task][training_set].append(rep_data)
                else:
                    try:
                        ref_acc_by_cluster[task][training_set].append(pd.read_excel(f"{data_dir}/train/replicate_models/train__{task}__{training_set}__{replicate}/allelic_imbalance/{file_name}", 
                                                                             sheet_name="roc_auc", index_col=0))
                    except:
                        ref_acc_by_cluster[task][training_set].append([])
    return ref_acc_by_cluster


def format_df_to_plot(accuracy_dict, cell_type_peak_set_mapping, tasks, training_sets):
    df_cols = ["Cell Type", "Cluster", "Training Task", "Cluster, Training Task", "Training Data", "Replicate", "Pearson R"]
    df_to_plot = pd.DataFrame([], columns=df_cols)
    
    for cell_type, clusters in cell_type_peak_set_mapping.items():
        for cluster in clusters:
            for task in tasks:
                for training_set in training_sets:
                    for replicate in range(1, 4):
                        try:
                            accuracy = accuracy_dict[task][training_set][replicate-1].loc[cluster, cell_type]
                        except:
                            accuracy = np.nan
                        df_to_plot = df_to_plot.append(pd.Series([cell_type, cluster, task, f"{cluster}, {task}", training_set, replicate, accuracy], index=df_cols), ignore_index=True)

    return df_to_plot


def read_multi_task_targets_predictions(model_dir, targets_dir, clusters, cell_types):
    preds_dict = {}

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

        preds_dict[cluster] = preds_df
    return preds_dict


def read_single_task_targets_predictions(model_dir, targets_dir, clusters, cell_types):
    preds_dict = {}

    for cluster in clusters:
        preds_df = []
        for cell_type in cell_types:
            preds = h5py.File(f"{model_dir.format(cell_type=cell_type)}/predict_beds/{cluster}/predict.h5", "r")
            preds_df.append(pd.DataFrame(np.nan_to_num(np.squeeze(preds["preds"][:,:,:])), columns=[f"{cell_type}_pred"]))
        preds_df = pd.concat(preds_df, axis=1)
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

        preds_dict[cluster] = preds_df
    return preds_dict


def read_biological_replicate_data(targets_dir, clusters, cell_types, replicates, chrs=["chr4", "chr5"]):
    preds_dict = {}
    for cluster in clusters:
        preds_df = pd.read_csv(f"{targets_dir}/{cluster}/predict_regions.bed", sep="\t", names=["chrom", "start", "end", "name"])

        for cell_type in cell_types:    
            for replicate in replicates:
                if os.path.exists(f"{targets_dir}/{cluster}/{cell_type}_{replicate}_target_signal.out"):
                    cell_type_targets = pd.read_csv(f"{targets_dir}/{cluster}/{cell_type}_{replicate}_target_signal.out", sep="\t", 
                                                   index_col=0, names=["size", "covered", "sum", "mean0", "mean"])
                    cell_type_targets = cell_type_targets.loc[preds_df["name"].values]["sum"].values
                    preds_df[f"{cell_type}_{replicate}_target"] = cell_type_targets
        
        preds_df = preds_df[preds_df["chrom"].isin(chrs)]
        preds_df = preds_df[(preds_df == 0).sum(axis=1) == 0]
        preds_dict[cluster] = preds_df
    return preds_dict


def compute_cross_cell_type_correlations(model_preds_dict, cell_type_peak_set_mapping, ubiquitous_cluster):
    correlations = pd.DataFrame([])
    for i, (task, preds_dict) in enumerate(model_preds_dict.items()):
        for cell_type, peak_sets in cell_type_peak_set_mapping.items():
            for peak_set in peak_sets:
                for cell_type_2 in cell_type_peak_set_mapping.keys():
                    if cell_type != cell_type_2:
                       
                        corr = pearsonr(np.log2(preds_dict[peak_set][f"{cell_type}_pred"]+1),
                                        np.log2(preds_dict[peak_set][f"{cell_type_2}_pred"]+1))[0]
                        correlations = correlations.append({"Training task": task,
                                                            "Cluster": "Ubiquitous" 
                                                            if peak_set == ubiquitous_cluster 
                                                            else "Cell-type specific",
                                                            "Cell Type": cell_type,
                                                            "Cell Type 2": cell_type_2,
                                                            "log-log Pearson R": corr}, ignore_index=True)
                        
                        if i == 0:
                            corr = pearsonr(np.log2(preds_dict[peak_set][f"{cell_type}_target"]+1),
                                            np.log2(preds_dict[peak_set][f"{cell_type_2}_target"]+1))[0]
                            correlations = correlations.append({"Training task": "Measurement",
                                                                "Cluster": "Ubiquitous" 
                                                                if peak_set == ubiquitous_cluster 
                                                                else "Cell-type specific",
                                                                "Cell Type": cell_type,
                                                                "Cell Type 2": cell_type_2,
                                                                "log-log Pearson R": corr}, ignore_index=True)
    return correlations 

def compute_biological_replicate_correlations(biological_replicate_data, cell_type_peak_set_mapping, replicates, ubiquitous_cluster):
    correlations = pd.DataFrame([])

    for cell_type, peak_sets in cell_type_peak_set_mapping.items():
        for peak_set in peak_sets:
            for cell_type_2 in cell_type_peak_set_mapping.keys():
                if cell_type != cell_type_2:
                    task = "Cross-individual measurement correlation (different cell types)"
                else:
                    task = "Cross-individual measurement correlation (same cell type)"
                for rep_1 in replicates:
                    for rep_2 in replicates:
                        if rep_1 != rep_2:
                            try:

                                corr = pearsonr(np.log2(biological_replicate_data[peak_set][f"{cell_type}_{rep_1}_target"]+1),
                                                np.log2(biological_replicate_data[peak_set][f"{cell_type_2}_{rep_2}_target"]+1))[0]
                                correlations = correlations.append({"Training task": task,
                                                                    "Cluster": "Ubiquitous" 
                                                                    if peak_set == ubiquitous_cluster 
                                                                    else "Cell-type specific",
                                                                    "Cell Type": cell_type,
                                                                    "Cell Type 2": cell_type_2,
                                                                    "log-log Pearson R": corr}, ignore_index=True)
                            except: 
                                continue
    return correlations 

def compute_false_positive_false_negative_error_rates(all_sequences, cell_type_peak_set_mapping, ubiquitous_cluster):
    error_df = pd.DataFrame([])

    for i, cell_type in enumerate(cell_type_peak_set_mapping.keys()):
        for j, tasks in enumerate(["Multi-task", "Single-task"]):
            # Total false neg rate
            fnr = len(all_sequences[(all_sequences[f"{cell_type} Peaks"]) &
                                    (~all_sequences[f"{tasks} {cell_type} predicted Peaks"])])/len(all_sequences[all_sequences[f"{cell_type} Peaks"]])

            # Total false pos rate
            fpr = len(all_sequences[(~all_sequences[f"{cell_type} Peaks"]) &
                                    (all_sequences[f"{tasks} {cell_type} predicted Peaks"])])/len(all_sequences[~all_sequences[f"{cell_type} Peaks"]])
            error_df = error_df.append({"Cluster": "All sequences",
                                        "Cell Type": cell_type,
                                        "Tasks": tasks,
                                        "False negative rate": fnr,
                                        "False positive rate": fpr}, ignore_index=True)

            # Peaks false negative rate
            fnr = len(all_sequences[(all_sequences[f"All peaks"] > 0) &
                                    (all_sequences[f"{cell_type} Peaks"]) &
                                    (~all_sequences[f"{tasks} {cell_type} predicted Peaks"])])/len(all_sequences[(all_sequences[f"All peaks"] > 0) &
                                                                                                                (all_sequences[f"{cell_type} Peaks"])])
            
            
            # Peaks false positive rate
            fpr = len(all_sequences[(all_sequences[f"All peaks"] > 0) &
                                    (~all_sequences[f"{cell_type} Peaks"]) &
                                    (all_sequences[f"{tasks} {cell_type} predicted Peaks"])])/len(all_sequences[(all_sequences[f"All peaks"] > 0) &
                                                                                                                (~all_sequences[f"{cell_type} Peaks"])])
            error_df = error_df.append({"Cluster": "All peaks",
                                        "Cell Type": cell_type,
                                        "Tasks": tasks,
                                        "False negative rate": fnr,
                                        "False positive rate": fpr}, ignore_index=True)

            for cluster in cell_type_peak_set_mapping[cell_type]:
                cluster_name = "Ubiquitous" if cluster == ubiquitous_cluster else "Cell type specific"
                fnr = len(all_sequences[(all_sequences[cluster] > 0) &
                          (~all_sequences[f"{tasks} {cell_type} predicted Peaks"])])/len(all_sequences[all_sequences[cluster] > 0])
                error_df = error_df.append({"Cluster": cluster,
                                        "Cell Type": cell_type,
                                        "Tasks": f"{tasks}, {cluster_name} peaks",
                                        "False negative rate": fnr,
                                        "False positive rate": np.nan}, ignore_index=True)
    return error_df
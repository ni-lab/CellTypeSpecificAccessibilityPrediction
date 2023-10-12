# CellTypeSpecificAccessibilityPrediction

This repository contains information and scripts to train and benchmark genomic deep learning models on their performance in cell-type specific regulatory regions. It makes use of pre-trained [Enformer](https://github.com/deepmind/deepmind-research/tree/master/enformer) and [Sei](https://github.com/FunctionLab/sei-manuscript) models, which can be downloaded from the linked sources. It also makes use of scripts from the [Basenji](https://github.com/calico/basenji/tree/master) repository for training additional models to benchmark different training decisions.

## Overview

We have organized the repository by figure. Within the `generate_figures/` directory, there is a notebook per figure to reproduce the analyses. Additional preprocessing and model training scripts are included in the relevant directories. 

This code has been tested on Python 3.7, and also makes use of Tensorflow 2.1. Please set up a conda environment and install the packages listed in the `requirements.txt` file.

## Data and models

Several different datasets and pre-trained models are used in these analyses. The following instructions can be used to download the relevant resources.

Much of the processed data and resources used in this repository can be found in the resources directory. Additional resources will be made available for download soon. This includes the model parameters and model weights of tissue-specific models, and various data used for evaluating performance in different peak regions. Details on how to download additional resources, such pre-trained models and large datasets hosted elsewhere, are described below.

#### Enformer

* The pre-trained Enformer model can be downloaded from TFhub ([link](https://tfhub.dev/deepmind/enformer/1))
* Enformer training, validation and test data can be downloaded from Google Cloud ([link](https://console.cloud.google.com/storage/browser/basenji_barnyard/data)). Note: This data is ~320 GB and is in a requester pays bucket.
* Pre-computed variant effect predictions for all frequent variants in the 1000 genomes cohort (MAF>0.05% in any population) can be downloaded from Google Cloud ([link](https://console.cloud.google.com/storage/browser/dm-enformer/variant-scores/1000-genomes/enformer;tab=objects?prefix=&forceOnObjectsSortingFiltering=false)). Note: this data is ~100GB.

#### Sei

* The pre-trained Sei model and relevant resources can be downloaded from Zenodo ([model](https://zenodo.org/record/4906997), [resources](https://zenodo.org/record/4906962))
* Sei test data and predictions can be downloaded from S3 using the folllowing instructions. Note: This data is 186 GB compressed and 612 GB decompressed.
  ```
  wget https://sei-files.s3.amazonaws.com/performance_curves.tar.gz
  tar -xzvf performance_curves.tar.gz
  ```

#### Calderon et al. data ([reference](https://www.nature.com/articles/s41588-019-0505-9))

* Bigwig files used to train models can be S3 using the `download_bigwigs.sh` script ([link](https://github.com/ni-lab/CellTypeSpecificAccessibilityPrediction/blob/main/scripts/tissue_specific_models/preprocess_calderon_data/download_bigwigs.sh))
* Cell type specific peaks and allelic imbalance data can be found in [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0505-9/MediaObjects/41588_2019_505_MOESM3_ESM.xlsx) of [Calderon et al. (2019)](https://www.nature.com/articles/s41588-019-0505-9). Cell type specific peaks are found in the sheet `lineage_groups` and allelic imbalance data is found in the sheet `significant_ASCs`.

#### Loeb et al. data 

* These data are available upon request and will be made public upon publication of Loeb, et al. Variants in tubule epithelial regulatory elements mediate most heritable differences in human kidney function. (Submitted).

#### Additional benchmark datasets

* GTeX SuSie fine-mapped eQTL data from [Wang et al. (2021)](https://www.nature.com/articles/s41592-021-01252-x#ref-CR22) and [Avsec et al. (2021)](https://www.nature.com/articles/s41592-021-01252-x) can be downloaded from Google Cloud ([link](https://console.cloud.google.com/storage/browser/dm-enformer/data/gtex_fine))
* UK Biobank GWAS summary statistics can be downloaded from the Neale lab server using the `download_gwas_sumstats.sh` script ([link](https://github.com/ni-lab/CellTypeSpecificAccessibilityPrediction/blob/main/scripts/enformer/ldsc/download_gwas_sumstats.sh))

## Analysis

Within the `scripts/` directory, a README within each subfolder describes how to perform the relevant analysis.

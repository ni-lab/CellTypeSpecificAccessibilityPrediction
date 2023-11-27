The scripts in this directory make use of a [forked version](https://github.com/pkathail/basenji/tree/snakemake/bin) of the basenji repository (pkathail/basenji:snakemake) to be able to train and evaluate models using a Snakemake pipeline. Once you have installed the forked version of the repository, the following scripts can be used to process datasets, train and evaluate models:
* `multi_task_preprocess_data.sh`: process datasets for all of the multi-task models
* `single_task_preprocess_data.sh`: process datasets for all of the single-task models
* `multi_task_Snakefile`: snakemake pipeline to train and evaluate the multi-task models
* `single_task_Snakefile`: snakemake pipeline to train and evaluate the single-task models

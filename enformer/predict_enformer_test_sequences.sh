#!/bin/bash
# Job name:
#SBATCH --job-name=enformer
#
# Account:
#SBATCH --account=co_nilah
#
# Partition:
#SBATCH --partition=savio3_gpu
#
# Quality of Service:
#SBATCH --qos=savio_lowprio
#
# Number of nodes:
#SBATCH --nodes=1
#
# Number of tasks (one for each GPU desired for use case) (example):
#SBATCH --ntasks=1
#
# Processors per task (please always specify the total number of processors twice the number of GPUs):
#SBATCH --cpus-per-task=4
#
#Number of GPUs, this can be in the format of "gpu:[1-4]", or "gpu:K80:[1-4] with the type included
#SBATCH --gres=gpu:1
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.kathail@berkeley.edu
#
## Command(s) to run:
module unload python
module unload cuda
module load ml/tensorflow/2.1.0-py37
export BASENJIDIR=/global/home/users/poojakathail/basenji
export PATH=$BASENJIDIR/bin:$PATH
export PYTHONPATH=$BASENJIDIR/bin:$PYTHONPATH
export TF_FORCE_GPU_ALLOW_GROWTH='true'

python /global/home/users/poojakathail/basenji/bin/sonnet_predict_bed.py --rc -a 2 -o /global/scratch/users/poojakathail/enformer/test -t /clusterfs/nilah/richard/basenji2/data_pretrained_basenji2/human/targets.txt -g /clusterfs/nilah/pooja/genomes/human.hg38.genome -f /clusterfs/nilah/pooja/genomes/hg38.ml.fa /global/scratch/users/poojakathail/enformer/saved_model/enformer /global/scratch/users/poojakathail/enformer/test_sequences_196608bp.bed

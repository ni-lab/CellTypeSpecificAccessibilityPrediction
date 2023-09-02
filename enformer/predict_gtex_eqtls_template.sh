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
#SBATCH --time=24:00:00
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.kathail@berkeley.edu
#
## Command(s) to run:
module unload python
module unload cuda
module load ml/tensorflow/2.5.0-py37
export BASENJIDIR=/global/home/users/poojakathail/basenji
export PATH=$BASENJIDIR/bin:$PATH
export PYTHONPATH=$BASENJIDIR/bin:$PYTHONPATH
export TF_FORCE_GPU_ALLOW_GROWTH='true'

# make predictions for variants in vcf files


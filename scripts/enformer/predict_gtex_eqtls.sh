#!/bin/bash

for tissue in `ls -1 /global/scratch/users/poojakathail/enformer/gtex_fine/susie | cut -d . -f1`
do
    echo $tissue
    cp predict_gtex_eqtls_template.sh enformer_pred_eqtl_vcf_${tissue}.sh
    for var_set in pos neg;
    do
        echo -e "python /global/home/users/poojakathail/basenji/bin/sonnet_sad.py --rc --shifts -1,0,1 -b 2 --stats SAD,REF,ALT -o /global/scratch/users/poojakathail/enformer/gtex_fine/preds/${tissue}_${var_set}_REF_ALT -t /clusterfs/nilah/richard/basenji2/data_pretrained_basenji2/human/targets.txt -f /clusterfs/nilah/pooja/genomes/hg38.ml.fa /global/scratch/users/poojakathail/enformer/saved_model/enformer /global/scratch/users/poojakathail/enformer/gtex_fine/vcf/${tissue}_${var_set}.vcf" >> enformer_pred_eqtl_vcf_${tissue}.sh
    done
    sbatch enformer_pred_eqtl_vcf_${tissue}.sh
done

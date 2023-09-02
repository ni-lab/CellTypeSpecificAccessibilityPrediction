#!/bin/bash

for tissue in `ls -1 /global/scratch/users/poojakathail/enformer/gtex_fine/susie | cut -d . -f1`
do
    echo $tissue
    cp predict_gtex_eqtls_template.sh sei_pred_eqtl_vcf_${tissue}.sh
    for var_set in pos neg;
    do
        echo -e "python /global/scratch/users/poojakathail/sei/sei-framework/1_variant_effect_prediction.py /global/scratch/users/poojakathail/enformer/gtex_fine/vcf/${tissue}_${var_set}.vcf /global/scratch/users/poojakathail/enformer/gtex_fine/preds/sei_${tissue}_${var_set} --genome=hg38 --cuda" >> sei_pred_eqtl_vcf_${tissue}.sh
    done
    sbatch sei_pred_eqtl_vcf_${tissue}.sh
done

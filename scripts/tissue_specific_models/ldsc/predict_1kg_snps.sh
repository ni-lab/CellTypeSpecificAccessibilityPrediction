#!/bin/bash

script_dir="/clusterfs/nilah/pooja/batch_scripts/sad"

dataset="immune_atlas"
for chr in {1..22}
do  
     tasks="multitask"
     mkdir /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__all_sequences__1
     cp ${script_dir}/immune_atlas_1kg_snps_sad_template.sh ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
     echo "python /global/home/users/poojakathail/basenji/bin/basenji_sad.py --stats SAD,REF,ALT -f /clusterfs/nilah/pooja/genomes/hg19.ml.fa -o /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__all_sequences__1/chr${chr} --rc --shifts 1,0,-1 -t /clusterfs/nilah/pooja/${dataset}/targets.txt /clusterfs/nilah/pooja/${dataset}/models/params_sc_kidney_regression_${tasks}.json /clusterfs/nilah/pooja/${dataset}/train/replicate_models/train__${tasks}__all_sequences__1/model_best.h5 /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink/1000G.EUR.QC.val.${chr}.vcf" >> ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
     sbatch ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
    
     tasks="single_task"
     for cell_type in DCmye NKim TCD8EM;
     do
         mkdir /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__${cell_type}__all_sequences__1
         cp ${script_dir}/immune_atlas_1kg_snps_sad_template.sh ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${cell_type}_${chr}.sh
         echo "python /global/home/users/poojakathail/basenji/bin/basenji_sad.py --stats SAD,REF,ALT -f /clusterfs/nilah/pooja/genomes/hg19.ml.fa -o /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__${cell_type}__all_sequences__1/chr${chr} --rc --shifts 1,0,-1 -t /clusterfs/nilah/pooja/${dataset}/single_task_targets/${cell_type}_targets.txt /clusterfs/nilah/pooja/${dataset}/models/params_sc_kidney_regression_${tasks}.json /clusterfs/nilah/pooja/${dataset}/train/replicate_models/train__${tasks}__${cell_type}__all_sequences__1/model_best.h5 /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink/1000G.EUR.QC.val.${chr}.vcf" >> ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${cell_type}_${chr}.sh
         sbatch ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${cell_type}_${chr}.sh
     done
    
    tasks="multitask_8x_params_same_layers"
    mkdir /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__all_sequences__1
    cp ${script_dir}/immune_atlas_1kg_snps_sad_template.sh ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
    echo "python /global/home/users/poojakathail/basenji/bin/basenji_sad.py --stats SAD,REF,ALT -f /clusterfs/nilah/pooja/genomes/hg19.ml.fa -o /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__all_sequences__1/chr${chr} --rc --shifts 1,0,-1 -t /clusterfs/nilah/pooja/${dataset}/targets.txt /clusterfs/nilah/pooja/${dataset}/models/params_sc_kidney_regression_${tasks}.json /clusterfs/nilah/pooja/${dataset}/train/replicate_models/train__${tasks}__all_sequences__1/model_best.h5 /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink/1000G.EUR.QC.val.${chr}.vcf" >> ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
    sbatch ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
done

dataset="kidney_data"
for chr in {1..22}
do
     tasks="multitask"
     mkdir /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__all_sequences__1
     cp ${script_dir}/immune_atlas_1kg_snps_sad_template.sh ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
     echo "python /global/home/users/poojakathail/basenji/bin/basenji_sad.py --stats SAD,REF,ALT -f /clusterfs/nilah/pooja/genomes/hg38.ml.fa -o /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__all_sequences__1/chr${chr} --rc --shifts 1,0,-1 -t /clusterfs/nilah/pooja/${dataset}/targets.txt /clusterfs/nilah/pooja/${dataset}/models/params_sc_kidney_regression_${tasks}.json /clusterfs/nilah/pooja/${dataset}/train/replicate_models/train__${tasks}__all_sequences__1/model_best.h5 /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink/1000G.EUR.QC.hg38.val.${chr}.vcf" >> ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
     sbatch ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
    
     tasks="single_task"
     for cell_type in PT LOH DT Str;
     do
         mkdir /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__${cell_type}__all_sequences__1
         cp ${script_dir}/immune_atlas_1kg_snps_sad_template.sh ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${cell_type}_${chr}.sh
         echo "python /global/home/users/poojakathail/basenji/bin/basenji_sad.py --stats SAD,REF,ALT -f /clusterfs/nilah/pooja/genomes/hg38.ml.fa -o /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__${cell_type}__all_sequences__1/chr${chr} --rc --shifts 1,0,-1 -t /clusterfs/nilah/pooja/${dataset}/single_task_targets/${cell_type}_targets.txt /clusterfs/nilah/pooja/${dataset}/models/params_sc_kidney_regression_${tasks}.json /clusterfs/nilah/pooja/${dataset}/train/replicate_models/train__${tasks}__${cell_type}__all_sequences__1/model_best.h5 /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink/1000G.EUR.QC.hg38.val.${chr}.vcf" >> ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${cell_type}_${chr}.sh
         sbatch ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${cell_type}_${chr}.sh
     done
    
    tasks="multitask_8x_params_same_layers"
    mkdir /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__all_sequences__1
    cp ${script_dir}/immune_atlas_1kg_snps_sad_template.sh ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
    echo "python /global/home/users/poojakathail/basenji/bin/basenji_sad.py --stats SAD,REF,ALT -f /clusterfs/nilah/pooja/genomes/hg38.ml.fa -o /clusterfs/nilah/pooja/${dataset}/ldsc/sad/train__${tasks}__all_sequences__1/chr${chr} --rc --shifts 1,0,-1 -t /clusterfs/nilah/pooja/${dataset}/targets.txt /clusterfs/nilah/pooja/${dataset}/models/params_sc_kidney_regression_${tasks}.json /clusterfs/nilah/pooja/${dataset}/train/replicate_models/train__${tasks}__all_sequences__1/model_best.h5 /clusterfs/nilah/pooja/immune_atlas/ldsc/geno/1000G_EUR_Phase3_plink/1000G.EUR.QC.hg38.val.${chr}.vcf" >> ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
    sbatch ${script_dir}/predict_1kg_snps_sad_${dataset}_${tasks}_${chr}.sh
done

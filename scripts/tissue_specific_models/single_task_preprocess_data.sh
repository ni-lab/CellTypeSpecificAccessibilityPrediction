#!/bin/bash

# preprocess datasets for all single-task models

outdir_prefix=$1
cell_types=$2
cell_type_list=$(echo $cell_types | tr "," "\n")
data_1_dir=$3
data_2_dir=$4
data_3_dir=$5
data_4_dir=$6
data_5_dir=$7
all_peaks=$8
non_ubiquitous_peaks=$9
targets_dir=$10
genome=$11

if [ $genome -eq "hg38" ]
then
	blocklist_file="/clusterfs/nilah/pooja/genomes/hg38.blacklist.rep.bed"
	gaps_file="/clusterfs/nilah/pooja/genomes/hg38_gaps.bed"
	genome_fasta="/clusterfs/nilah/pooja/genomes/hg38.ml.fa"
	genome_sizes="/clusterfs/nilah/pooja/genomes/human.hg38.genome"
else
	blocklist_file="/clusterfs/nilah/pooja/genomes/wgEncodeHg19ConsensusSignalArtifactRegions.bed"
	gaps_file="/clusterfs/nilah/pooja/genomes/hg19_gaps.bed"
	genome_fasta="/clusterfs/nilah/pooja/genomes/hg19.ml.fa"
	genome_sizes="/clusterfs/nilah/pooja/genomes/human.hg19.genome"
fi

for cell_type in $cell_type_list;
do
	echo ${cell_type}
	
	###
	# 1. All sequences
	###
	python /global/home/users/poojakathail/basenji/bin/basenji_data.py -b ${blocklist_file} -g ${gaps_file} -p 20 -r 4096 -w 192 -l 1344 -v "chr7,chr14,chr15" -t "chr4,chr5" --stride 192 --stride_test 192 --crop 576 -o ${outdir_prefix}/${cell_type}/${data_1_dir} --local ${genome_fasta} ${targets_dir}/${cell_type}_targets.txt

	###
	# 2. All sequences overlapping a peak
	###
	python /global/home/users/poojakathail/basenji/bin/basenji_data.py --limit ${all_peaks} -b ${blocklist_file} -g ${gaps_file} -p 20 -r 4096 -w 192 -l 1344 -v "chr7,chr14,chr15" -t "chr4,chr5" --stride 192 --stride_test 192 --crop 576 -o ${outdir_prefix}/${cell_type}/${data_2_dir} --local ${genome_fasta} ${targets_dir}/${cell_type}_targets.txt

	###
	# 3. All sequences overlapping a non-ubiquitous peak
	###
	python /global/home/users/poojakathail/basenji/bin/basenji_data.py --limit ${non_ubiquitous_peaks} -b ${blocklist_file} -g ${gaps_file} -p 20 -r 4096 -w 192 -l 1344 -v "chr7,chr14,chr15" -t "chr4,chr5" --stride 192 --stride_test 192 --crop 576 -o ${outdir_prefix}/${cell_type}/${data_3_dir} --local ${genome_fasta} ${targets_dir}/${cell_type}_targets.txt

	# extract sequences that don't overlap peaks
	bedtools intersect -a ${outdir_prefix}/${cell_type}/${data_1_dir}/sequences.bed -b ${all_peaks} -v > ${outdir_prefix}/${cell_type}/non_peaks.bed

	# annotate peak and non-peak sequences with GC content
	bedtools nuc -fi ${genome_fasta} -bed ${outdir_prefix}/${cell_type}/${data_2_dir}/sequences.bed > ${outdir_prefix}/${cell_type}/${data_2_dir}/sequences_nuc.bed
	bedtools nuc -fi ${genome_fasta} -bed ${outdir_prefix}/${cell_type}/non_peaks.bed > ${outdir_prefix}/${cell_type}/non_peaks_nuc.bed

	###
	# 4. 1:1 peaks:non-peaks
	###

	# generate balanced training sequences
	python generate_balanced_training_sequences.py --peaks ${outdir_prefix}/${cell_type}/${data_2_dir}/sequences_nuc.bed --non_peaks ${outdir_prefix}/${cell_type}/non_peaks_nuc.bed --out ${outdir_prefix}/${cell_type}/${data_4_dir}

	python /global/home/users/poojakathail/basenji/bin/basenji_data.py --restart -b ${blocklist_file} -g ${gaps_file} -p 20 -r 4096 -w 192 -l 1344 -v "chr7,chr14,chr15" -t "chr4,chr5" --stride 192 --stride_test 192 --crop 576 -o ${outdir_prefix}/${cell_type}/${data_4_dir} --local ${genome_fasta} ${targets_dir}/${cell_type}_targets.txt

	###
	# 5. 1:1 peaks:non-peaks (GC matched)
	###

	# generate balanced training sequences
	python generate_balanced_training_sequences.py --peaks ${outdir_prefix}/${cell_type}/${data_2_dir}/sequences_nuc.bed --non_peaks ${outdir_prefix}/${cell_type}/non_peaks_nuc.bed --out ${outdir_prefix}/${cell_type}/${data_5_dir} --gc

	python /global/home/users/poojakathail/basenji/bin/basenji_data.py --restart -b ${blocklist_file} -g ${gaps_file} -p 20 -r 4096 -w 192 -l 1344 -v "chr7,chr14,chr15" -t "chr4,chr5" --stride 192 --stride_test 192 --crop 576 -o ${outdir_prefix}/${cell_type}/${data_5_dir} --local ${genome_fasta} ${targets_dir}/${cell_type}_targets.txt
done



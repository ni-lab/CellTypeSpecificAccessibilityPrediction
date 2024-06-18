#!/bin/bash

# only needs to be run for Loeb et al. peaks which are in hg38

ubiquitous_peaks_prefix=$1
cell_type_specific_peaks_prefix=$2

/clusterfs/nilah/pooja/software/liftOver_linux ${ubiquitous_peaks_prefix}.bed /clusterfs/nilah/pooja/software/hg38ToHg19.over.chain ${ubiquitous_peaks_prefix}.hg19.bed ${ubiquitous_peaks_prefix}.unmap

/clusterfs/nilah/pooja/software/liftOver_linux ${cell_type_specific_peaks_prefix}.bed /clusterfs/nilah/pooja/software/hg38ToHg19.over.chain ${cell_type_specific_peaks_prefix}.hg19.bed ${cell_type_specific_peaks_prefix}.unmap

/clusterfs/nilah/pooja/software/liftOver_linux cluster11_PT.bed /clusterfs/nilah/pooja/software/hg38ToHg19.over.chain cluster11_PT.hg19.bed cluster11_PT.unmap
/clusterfs/nilah/pooja/software/liftOver_linux cluster4_DistalNephron.bed /clusterfs/nilah/pooja/software/hg38ToHg19.over.chain cluster4_DistalNephron.hg19.bed cluster4_DistalNephron.unmap
/clusterfs/nilah/pooja/software/liftOver_linux cluster9_Stroma.bed /clusterfs/nilah/pooja/software/hg38ToHg19.over.chain cluster9_Stroma.hg19.bed cluster9_Stroma.unmap




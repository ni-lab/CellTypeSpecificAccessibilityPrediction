#!/bin/bash

# ls -1 *.bw | cut -d '_' -f1 | uniq > cell_types.txt
while read p; do
  /clusterfs/nilah/pooja/software/bigWigMerge ${p}_*.bw ../merged/${p}_merged.bedGraph
  sort -k1,1 -k2,2n ../merged/${p}_merged.bedGraph > ../merged/${p}_merged.sorted.bedGraph
  /clusterfs/nilah/pooja/software/bedGraphToBigWig ../merged/${p}_merged.sorted.bedGraph /clusterfs/nilah/pooja/genomes/human.hg19.genome ../merged/${p}_merged.bw
  rm ../merged/*.bedGraph
done <cell_types.txt

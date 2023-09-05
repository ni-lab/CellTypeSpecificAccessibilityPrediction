#!/bin/bash

wget https://s3.amazonaws.com/muellerf/data/trackhubs/immune_atlas/hg19/trackDb.txt
grep .bw trackDb.txt | cut -d ' ' -f2 > bigwig_files.txt
while read p; do
  wget https://s3.amazonaws.com/muellerf/data/trackhubs/immune_atlas/hg19/${p}
done <bigwig_files.txt


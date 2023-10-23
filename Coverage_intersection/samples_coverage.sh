#!/bin/bash

# Generate the case bed file where 90% samples have coverage >= 10 

outputfile="summary.tsv.gz"

zcat "$outputfile" | awk '(NR > 1 && $7 >= 0.9) {print $1,$2 - 1,$2}' | sort -k1,1 -k2,2n | sed 's/ /\t/g' | bedtools merge -i stdin | gzip > sample.coverage10.bed.gz

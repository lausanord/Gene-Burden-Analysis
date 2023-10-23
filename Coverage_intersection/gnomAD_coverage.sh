#!/bin/bash

# Restrict gnomAD coverage data to the region where >= 90% of samples have coverage of depth >=10.

gnomADCoverage="gnomad.exomes.r2.1.1.sites.vcf.gz"
bedFile="gnomad.coverage10.bez.gz"

zcat ${gnomADCoverage} | awk '(NR > 1 && $7 >= 0.9) {print $1,$2 - 1,$2}' \
  | sort -k1,1 -k2,2n | sed 's/ /\t/g' | bedtools merge -i stdin | gzip > ${bedFile}

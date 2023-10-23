#!/bin/bash

# 1: Compress all VCF files
for file in *.vcf; do
    bgzip "$file"
done

# 2: Index compressed files
for file in *.vcf.gz; do
    bcftools index "$file"
done

# 3: Creates a list of compressed VCF files
find . -name "*.vcf.gz" > samples.txt

# 4: Joins VCF files into one
bcftools merge -l samples.txt -O v -o merge.vcf.gz

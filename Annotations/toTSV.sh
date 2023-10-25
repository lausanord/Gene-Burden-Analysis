#!/bin/bash

# Take a VCF file to convert it to a TSV file.

## Define input and output filename as variables.
input_vcf="samples.vcf"
output_file="output.tsv"

## Process input file and add results to output file.
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\n' -d -A tab $input_vcf

echo -e "CHROM\tPOS\tREF\tALT\t$(bcftools +split-vep -l $input_vcf | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > $output_file

bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\n' -d -A tab $input_vcf >> $output_file

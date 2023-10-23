#!/bin/bash
set -eu -o pipefail

main() {
  ## split multiallelic variants into biallelic variants, left normalization,
  ## and QC on DP, GQ, and VAF for heterozygous variants 

  ## input
  # vcfFile=$1
  # outputPrefix=$2
  # refFASTA=$3
  # BCFTOOLS=$4

  ## output
  ## ${outputPrefix}.biallelic.leftnorm.ABCheck.vcf.gz: vcf gz file
  ## ${outputPrefix}.biallelic.leftnorm.ABCheck.vcf.gz.tbi: index file
  ## ${outputPrefix}.biallelic.leftnorm.filtered.txt: the failed variants from
  ## the Filter column


  vcfFile=$1
  outputPrefix=$2
  refFASTA=$3
  BCFTOOLS=$4
  
  mkdir -p $(dirname ${outputPrefix})

  ## biallelic, normalize and 
  # AB filtering: for het, VAF >= 0.2 & VAF <= 0.8
  ${BCFTOOLS} norm -m-both -f ${refFASTA} -Ou ${vcfFile} | \
  ${BCFTOOLS} filter -S. -i \
    '(FMT/DP>=10 & FMT/GQ>=20) & (GT="hom" | (GT="het" & FMT/AD[:1] > 0 & (FMT/AD[:1] / (FMT/AD[:0] + FMT/AD[:1]) >= 0.2 - 1e-6) & (FMT/AD[:1] / (FMT/AD[:0] + FMT/AD[:1]) <= 0.8) ))' \
  -Oz -o \
  ${outputPrefix}.biallelic.leftnorm.ABCheck.vcf.gz 
  tabix -p vcf \
   ${outputPrefix}.biallelic.leftnorm.ABCheck.vcf.gz

  # store those variants that are filtered by VQSR
  vqsrFiltered=${outputPrefix}.biallelic.leftnorm.filtered.txt
  ${BCFTOOLS} query -f '%CHROM-%POS-%REF-%ALT\n' -i '%FILTER!="." & %FILTER!="PASS"' ${outputPrefix}.biallelic.leftnorm.ABCheck.vcf.gz > ${vqsrFiltered}
}

main "$@"
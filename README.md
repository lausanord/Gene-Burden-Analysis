# Gene-Burden-Analysis

#### We show the development of different tools to obtain the necessary files to carry out the CoCoRV (Consistent summary Count based Rare Variant burden test) framework, developed by Wenan Chen and Saima Sultana Tithi.

## 1. Samples preparation
Once we have our cases stored in VCF files, we proceed to aggregate the individual VCF files into a single comprehensive VCF file using the merge_VCF.sh script. 

As part of data quality control, we perform VCF normalization using the vcfQCAndNormalize.sh script. In some cases, VCF files may contain multiallelic variants, which can complicate subsequent analyses. To address this, we convert multiallelic variants into biallelic format, ensuring that the data is amenable to a wide range of genetic analyses.

## 2. Annotation
To annotate variants using ANNOVAR, we employ the annotate.sh script. This script is executed separately for cases and controls, facilitated by the annovar_cases.sh and annovar_controls.sh scripts, respectively.
To annotate variants using VEP, we use the VEPAnnotation.sh script.

For the consequence analysis, we use the script count_csq_one_sample.R to count how many variants there are of each type. To count the possibly pathogenic variants with each prediction algorithm we use toCSV.sh to obtain the table of the VCF and count_consequences.R.

The annotated VCF files, both cases and controls, can be converted to GDS format using the vcf2gds.R script.

## 3. Intersection coverage
To generate the intersection file of coverage of samples and controls in which at least 90% of the samples have coverage with depth 10, we convert bam to bw with the scripts bamToBW.sh and BW.sh, and extract the coverage summaries with bwToCoverageSummary.py and Coverag_summaries.sh. The summaries are merged with Summary_files.sh

The restricted coverages are extracted with samples_coverage.sh and gnomAD_coverage.sh and we do the intersection with intersection.sh.

## 4. Burden Test
We include the necessary R scripts for the designs with REVEL, CADD and AlphaMissense, in order to select the variants of interest in each design. 

In test.sh it is necessary to indicate the parameters and files to be included for the burden test.

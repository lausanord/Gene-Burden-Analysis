#!/bin/bash
VEP_CACHE="$HOME/vep_data/"
VEP_FASTA="$HOME/References/vep/VEPfasta/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
PLUGIN_DBS="$HOME/References/vep/VEPdbs"

singularity run /mnt/beegfs/home/mimarbor/Laura/PruebaVEPRevel/vep.sif vep \
-i $1 -o $1.VEP.vcf.gz \
--offline \
--fasta $VEP_FASTA \
--dir $VEP_CACHE \
--assembly GRCh37 \
--species homo_sapiens \
--cache \
--pick \
--vcf \
--no_stats \
--format vcf \
--plugin REVEL,file=$HOME/PruebaVEPRevel/new_tabbed_revel.tsv.gz,no_match=1 \
--plugin CADD,${PLUGIN_DBS}/InDels.tsv.gz,${PLUGIN_DBS}/whole_genome_SNVs.tsv.gz \
--plugin AlphaMissense,file=%HOME/AlphaMissense_hg19.tsv.gz \
--compress_output bgzip

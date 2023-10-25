# Load the vcfR package
library(vcfR)
# Load the stringr package
library(stringr)

# Path to the folder containing VCF files


# Initialize a dataframe to store the results
total_results = data.frame(Variant_Type = character(), Frequency = integer(), stringsAsFactors = FALSE)
total_variants = data.frame(Sample = character(), Variants = integer(), stringsAsFactors = FALSE)


  # Read the VCF file
  vcf_data = read.vcfR("samples.annotated.biallelic.leftnorm.ABCheck.noCHR.vcf.gz.VEP.vcf", verbose = FALSE)
  
  # Extract information from CSQ and keep only the variant types
  info_CSQ = extract_info_tidy(vcf_data, info_fields = "CSQ", info_types = TRUE, info_sep = ";")
  matches = str_match(info_CSQ$CSQ, "\\|([^|]+)\\|")
  info_CSQ$CSQ = matches[, 2]
  
  # Count the frequency of each unique value in the CSQ column
  frequency_elements = table(info_CSQ$CSQ)
  variant_type_count = as.data.frame(frequency_elements)
  variant_type_count$File = basename("samples.annotated.biallelic.leftnorm.ABCheck.noCHR.vcf.gz.VEP.vcf")  # Add a column for the file name
  
  # Rename the columns
  colnames(variant_type_count) = c("Type", "Count", "File")
  
  # Add the results to the total dataframe
  total_results = rbind(total_results, variant_type_count)
  
  # Count the number of variants
  variants_number = nrow(info_CSQ)
  total_variants = rbind(total_variants, data.frame(Sample = basename("samples.annotated.biallelic.leftnorm.ABCheck.noCHR.vcf.gz.VEP.vcf"), Variants = variants_number))


### Total variants
# Print the total variants for each file
print(total_variants)

# Print the sum of all total variants
total_variants_sum = sum(total_variants$Variants)
print(total_variants_sum)


### Variant type
# Print the different types of variants in each sample
print(total_results)

# Aggregate the counts across all files
final_results = aggregate(Count ~ Type, data = total_results, sum)

# Print the sum of all different types of variants
print(final_results)


### Aggregation of the variants in groups for each sample and total
## 5'-UTR
types_5UTR = c("5_prime_UTR_variant", "splice_region_variant&5_prime_UTR_variant")
count_5UTR = subset(total_results, Type %in% types_5UTR)
sum(count_5UTR$Count)

## 3'-UTR
types_3UTR = c("3_prime_UTR_variant","3_prime_UTR_variant&NMD_transcript_variant")
count_3UTR = subset(total_results, Type %in% types_3UTR)
sum(count_3UTR$Count)

## Exonic
# Nonsense
types_nonsense = c("stop_gained","stop_gained&frameshift_variant&splice_region_variant","stop_gained&frameshift_variant, stop_gained&splice_region_varian, stop_retained_variant", "start_lost", "stop_gained&splice_region_variant", "stop_lost")
count_nonsense = subset(total_results, Type %in% types_nonsense)
sum(count_nonsense$Count)

# INDELs
types_indel = c("frameshift_variant", "frameshift_variant&stop_lost",  "inframe_deletion", "inframe_insertion", "inframe_deletion&splice_region_variant", "inframe_insertion&splice_region_variant", "splice_acceptor_variant
", "splice_acceptor_variant&coding_sequence_variant", "	
splice_acceptor_variant&coding_sequence_variant&5_prime_UTR_variant&intron_variant", "splice_acceptor_variant&coding_sequence_variant&intron_variant", "splice_acceptor_variant&splice_polypyrimidine_tract_variant&intron_variant", "splice_donor_5th_base_variant&intron_variant", "splice_donor_region_variant&intron_variant" )
count_indels = subset(total_results, Type %in% types_indel)
sum(count_indels$Count)

# No coding
types_no_coding = c("non_coding_transcript_exon_variant", "intron_variant&non_coding_transcript_variant",  "splice_region_variant&non_coding_transcript_exon_variant", "splice_acceptor_variant&non_coding_transcript_variant", "splice_donor_5th_base_variant&intron_variant&non_coding_transcript_variant", "splice_donor_region_variant&intron_variant&non_coding_transcript_variant", "splice_region_variant&intron_variant&non_coding_transcript_variant", "splice_region_variant&non_coding_transcript_variant", "splice_region_variant&splice_polypyrimidine_tract_variant&intron_variant&non_coding_transcript_variant
")
count_no_coding = subset(total_results, Type %in% types_no_coding)
sum(count_no_coding$Count)

## Intronic
types_intronic = c("splice_donor_variant", "intron_variant", "splice_polypyrimidine_tract_variant&intron_variant", "splice_region_variant&intron_variant",  "splice_acceptor_variant&non_coding_transcript_variant", "splice_donor_variant&non_coding_transcript_variant", "splice_acceptor_variant&frameshift_variant", "intron_variant&NMD_transcript_variant", "splice_donor_variant&NMD_transcript_variant", "splice_donor_variant&splice_donor_5th_base_variant&coding_sequence_variant&intron_variant", "splice_donor_variant&splice_donor_region_variant&coding_sequence_variant&intron_variant", "splice_donor_variant&splice_donor_region_variant&intron_variant", "splice_polypyrimidine_tract_variant&intron_variant&non_coding_transcript_variant", "splice_region_variant&intron_variant", "splice_region_variant&splice_polypyrimidine_tract_variant&intron_variant")
count_intronic = subset(total_results, Type %in% types_intronic)
sum(count_intronic$Count)

## Intergenic
types_intergenic= c("intergenic_variant", "downstream_gene_variant", "upstream_gene_variant", "TF_binding_site_variant", "regulatory_region_variant", "TFBS_ablation&TF_binding_site_variant")
count_intergenic = subset(total_results, Type %in% types_intergenic)
sum(count_intergenic$Count)

## Intergenic
types_missense= c("missense")
count_missense = subset(total_results, Type %in% types_missense)
sum(count_missense$Count)

total_sum = sum(count_5UTR$Count, count_3UTR$Count, count_nonsense$Count, count_indels$Count, count_no_coding$Count, count_intronic$Count, count_intergenic$Count)
print(total_sum)









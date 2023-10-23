library(SeqArray)

main = function() {
  # convert vcf to gds 
  
  args = commandArgs(T)
  
  vcfFile = args[1]
  gdsFile = args[2]
  
  ncore = 1
  
  if (length(args) >= 3) {
    ncore = as.integer(args[3])
  }
  
  stopifnot(ncore >= 1)
  
  seqVCF2GDS(vcfFile, gdsFile, parallel = ncore, 
             ignore.chr.prefix="")
}

main()

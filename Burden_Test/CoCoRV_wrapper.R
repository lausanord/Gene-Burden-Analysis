suppressPackageStartupMessages(require("argparse"))

require(CoCoRV)

main = function() {
  ## association test using summary counts
  ## if stratified analysis is specified, the AC and AN will still be used for
  ## filtering and QC, then the counts corresponding to the stratified groups
  ## will be used in the CMH exact test
  
  parser = ArgumentParser()
  parser$add_argument("controlGDSFile", 
                      help = "VCF format file of controls")
  parser$add_argument("caseGDSFile", 
                      help = "VCF format file of cases")
  parser$add_argument("--controlAnnoGDSFile", metavar='',
                      action = "store", default = "",
                      help = "VCF format file of annotated variants in controls")
  parser$add_argument("--caseAnnoGDSFile", metavar='',
                      action = "store", default = "",
                      help = "VCF format file of annotated variants in cases")
  parser$add_argument("--nControl", metavar='',
                      action = "store", type = "integer", default = -1,
                      help = "number of controls, will be detected based on AN if not set")
  parser$add_argument("--nCase", metavar='',
                      action = "store", type = "integer", default = -1,
                      help = "number of cases, will be detected based on AN if not set")
  parser$add_argument("--sampleList", metavar='',
                      action = "store", default = "",
                      help = "a one column file listing sample IDs used in cases")
  parser$add_argument("--outputPrefix", metavar='',
                      action = "store", default = "out",
                      help = "the output file prefix")
  parser$add_argument("--AFMax", default = 1e-3, metavar='',
                      action = "store", type = "double",
                      help = "the maximum of the alternate alleles frequency")
  parser$add_argument("--AFUse", metavar='',
                      action = "store", default = "case_control_joint",
                      help = paste("how to apply the AF threshold, can be",
                                   "case_control_joint (default), control, case_control_separate"))
  parser$add_argument("--AFGroup", metavar='',
                      action = "store", default = "max",
                      help = paste("how to combine AFs from multiple groups, can be",
                                   "max: use the max of multiple group AFs (default)", 
                        ", pool: pool all groups together to calculate the AF"))
  parser$add_argument("--countUse", metavar='',
                      action = "store", default = "sample",
                      help = paste("how to use the count for association test",
              "sample (default)", 
    "allele_null (functional allele counts vs nonfunctional allele counts as in proxecat)"))
  parser$add_argument("--variantMissing", metavar='',
                      action = "store", default = 0.1, type = "double",
                      help = paste("the maximum missingness allowed for a variant"))
  parser$add_argument("--removeStar", 
                      action = "store_true", default = F,
                      help = paste("remove overlapping deletion variants with the star symbol.",
                                   "This is useful when using gnomAD because gnomAD does not have variants with",
                                   "star as alternate allele."))
  parser$add_argument("--groupColumn", metavar='',
                      action = "store", default = "Gene.refGene",
                      help = paste("the column specifying the column name grouping variants,",
                                   "such as the gene name"))
  parser$add_argument("--ACANConfig", metavar='',
                      action = "store", default = "",
                      help = paste("a file specifying the AC AN info ID and optional",
                                   "stratified ID for cases"))
  parser$add_argument("--caseGroup", metavar='',
                      action = "store", default = "",
                      help = paste("a file with sample ID in the first column and the",
                                   "stratified group info in the 2nd column"))
  parser$add_argument("--variantExcludeFile", metavar='',
                      action = "store", default = "",
                      help = paste("a one column file specifying the variants to be excluded"))
  parser$add_argument("--ignoreFilter", 
                      action = "store_true", default = F,
                      help = paste("if not specified, remove variants with filter not equal", 
    "to PASS in either case or control for a consistent filtering"))
  parser$add_argument("--variantGroupCustom", metavar='',
                      action = "store", default = "",
      help = paste("an optional R file with functions or a tab separated two column text file defining the variants of interest, see README for examples"))
  parser$add_argument("--extraParamJason", metavar='',
                      action = "store", default = "",
      help = paste("extraParamJason an optional JASON file to provide extra parameters to the custuomized file variantGroupCustom"))
  parser$add_argument("--variantGroup", metavar='',
                      action = "store", default = "pathogenic",
                      help = paste("a specified variant group to use for test or a self defined",
                                   "function to define the variants of interest, see the function CoCoRV for",
                                   "more details"))
  parser$add_argument("--annotationList", metavar='',
                      action = "store", default = "",
            help = paste("a one column file without the header listing all required annotations used in variant filtering. AC AN annotations from ACANConfig will be added automatically"))
  parser$add_argument("--minREVEL", metavar='',
                      action = "store", default = 0.65, type = "double",
                      help = paste("the minimum REVEL score for pathogenic missense variants"))
  parser$add_argument("--maxAFPopmax", metavar='',
                      action = "store", default = 1, type = "double",
                      help = paste("the maximal AF_popmax for filtering variants.",
                                   "for consistency, remove the AF_popmax from gnomAD raw vcf file and then",
                                   "annotate using another annotation tool, e.g., ANNOVAR"))
  parser$add_argument("--fullCaseGenotype", 
                      action = "store_true", default = F,
                      help = paste("If true (default), use the full genotypes to extract",
                                   "counts for cases, otherwise summary counts are used"))
  parser$add_argument("--bed", metavar='',
                      action = "store", default = "",
                      help = paste("a bed file to specify the regions considered, e.g., the",
                                   "intersection of the bed files from cases and controls,", 
                                   "where >= 90 percent of samples have coverage >= 10"))
  parser$add_argument("--overlapType", metavar='',
                      action = "store", default = "within",
                      help = paste("overlap type between the variant and the bed file:",
                                   "any, within (default)"))
  parser$add_argument("--checkHighLDInControl", 
                      action = "store_true", default = F,
                      help = paste("for each case double het, if high LD in control then exclude"))
  parser$add_argument("--pLDControl", metavar='',
                      action = "store", default = 0.05, type = "double",
                      help = paste("the pvalue threshold to detect high LD variants in control"))
  parser$add_argument("--ORThresholdControl", metavar='',
                      action = "store", default = 1, type = "double",
                      help = paste("the OR threshold to define high LD variants"))
  parser$add_argument("--highLDVariantFile", metavar='',
                      action = "store", default = "",
                      help = paste("precomputed high LD variant pairs in controls"))
  parser$add_argument("--ignoreEthnicityInLD", 
                      action = "store_true", default = T,
                      help = paste("whether to ignore ethnicity when loading the", 
                                   "precomputed LD pairs"))
  parser$add_argument("--gnomADVersion", metavar='',
                      action = "store", default = "v2exome",
                      help = paste("if gnomAD is used, specify the version:", 
                                   "v2exome, v2genome, v3genome"))
  parser$add_argument("--regionFile", metavar='',
                      action = "store", default = "",
                      help = paste("a bed file with names at the 4th column defining",  
                                   "the processing unit, can be used to define test regions or to", 
                                   "process relatively small regions to reduce memory"))
  parser$add_argument("--reference", metavar='',
                      action = "store", default = "GRCh37",
                      help = paste("the reference build, either GRCh37 or GRCh38,", "used for defining the pseudoautosomal region in X and Y, default GRCh37"))
  
  tryCatch({arguments = parser$parse_args()},
           error = function(e) {
             print(e)
             stop()
           })
  
  nControl = arguments$nControl
  nCase = arguments$nCase
  controlGDSFile = arguments$controlGDSFile
  caseGDSFile = arguments$caseGDSFile
  controlAnnoGDSFile = arguments$controlAnnoGDSFile
  caseAnnoGDSFile = arguments$caseAnnoGDSFile
  sampleList = arguments$sampleList
  outputPrefix = arguments$outputPrefix
  AFMax = arguments$AFMax
  AFUse = arguments$AFUse
  AFGroup = arguments$AFGroup
  countUse = arguments$countUse
  variantMissing = arguments$variantMissing
  removeStar = arguments$removeStar
  minREVEL = arguments$minREVEL
  maxAFPopmax = arguments$maxAFPopmax
  ACANConfig = arguments$ACANConfig
  caseGroup = arguments$caseGroup
  variantExcludeFile = arguments$variantExcludeFile
  ignoreFilter = arguments$ignoreFilter
  variantGroupCustom = arguments$variantGroupCustom
  extraParamJason = arguments$extraParamJason
  variantGroup = arguments$variantGroup
  groupColumn = arguments$groupColumn
  annotationList = arguments$annotationList
  bed = arguments$bed
  fullCaseGenotype = arguments$fullCaseGenotype
  overlapType = arguments$overlapType
  pLDControl = arguments$pLDControl
  ORThresholdControl = arguments$ORThresholdControl
  checkHighLDInControl = arguments$checkHighLDInControl
  highLDVariantFile = arguments$highLDVariantFile
  ignoreEthnicityInLD = arguments$ignoreEthnicityInLD
  gnomADVersion = arguments$gnomADVersion
  regionFile = arguments$regionFile
  
  argumentNames = names(arguments)
  cat("Input arguments:\n")
  for (i in 1:length(argumentNames)) {
    cat(paste0(argumentNames[i], "=", arguments[[i]]), "\n")
  }
  
  browser()
  resultPerRegion = CoCoRV(controlGDSFile = controlGDSFile,
                           caseGDSFile = caseGDSFile,
                           controlAnnoGDSFile = controlAnnoGDSFile,
                           caseAnnoGDSFile = caseAnnoGDSFile,
                           nControl = nControl,
                           nCase = nCase,
                           sampleList = sampleList,
                           outputPrefix = outputPrefix,
                           AFMax = AFMax,
                           AFUse = AFUse,
                           AFGroup = AFGroup,
                           countUse = countUse,
                           variantMissing = variantMissing,
                           removeStar = removeStar,
                           minREVEL = minREVEL,
                           maxAFPopmax = maxAFPopmax,
                           ACANConfig = ACANConfig,
                           caseGroup = caseGroup,
                           variantExcludeFile = variantExcludeFile,
                           ignoreFilter = ignoreFilter,
                           variantGroupCustom = variantGroupCustom,
                           extraParamJason = extraParamJason,
                           variantGroup = variantGroup,
                           groupColumn = groupColumn,
                           annotationList = annotationList,
                           bed = bed,
                           fullCaseGenotype = fullCaseGenotype,
                           overlapType = overlapType,
                           pLDControl = pLDControl,
                           ORThresholdControl = ORThresholdControl,
                           checkHighLDInControl = checkHighLDInControl,
                           gnomADVersion = gnomADVersion,
                           highLDVariantFile = highLDVariantFile,
                           ignoreEthnicityInLD = ignoreEthnicityInLD,
                           regionFile = regionFile)
  

  if (length(resultPerRegion) > 0) {
    file.create(paste0(outputPrefix, ".case.group"))
    file.create(paste0(outputPrefix, ".control.group"))
    caseGroup = file(paste0(outputPrefix, ".case.group"), open = "a")
    controlGroup = file(paste0(outputPrefix, ".control.group"), open = "a")
    for (i in 1:length(resultPerRegion)) {
      result = resultPerRegion[[i]]
      regionID = result$regionID
      
      if (regexpr("sample", countUse) != -1) {
        # output the group list files
        caseGroupLines = sapply(1:length(result$caseGeneList),
                                function(x) paste(result$caseGeneList[[x]], collapse=" "))
        controlGroupLines = sapply(1:length(result$controlGeneList),
                                   function(x) paste(result$controlGeneList[[x]], 
                                                     collapse=" "))
        if (!is.null(regionID)) {
          caseGroupLines = paste(regionID, caseGroupLines)
          controlGroupLines = paste(regionID, controlGroupLines)
          result$testResult = cbind(regionID, result$testResult)
        }
        writeLines(caseGroupLines, con = caseGroup)
        writeLines(controlGroupLines, con = controlGroup)
      }
      
      # output the test result
      if (i == 1) {
        write.table(result$testResult, 
                    file = paste0(outputPrefix, ".association.tsv"), 
                    row.names = F, col.names = T, sep = "\t", quote = F)
      } else {
        write.table(result$testResult, 
                    file = paste0(outputPrefix, ".association.tsv"), 
                    row.names = F, col.names = F, sep = "\t", quote = F, 
                    append = T)
      }
    }
    close(caseGroup)
    close(controlGroup)
  }
  
  cat("Analysis Done!\n")
}

main()

suppressPackageStartupMessages(require("argparse"))

library(CoCoRV)

main = function() {
  ## QQ plot for count tables and FDR calculation 
  
  parser = ArgumentParser()
  parser$add_argument("inputFile", 
          help = paste("The input file, with header, including at least one",
    "column of p values or 2x2xm counts, and one column of ID for each row"))
  parser$add_argument("outputPrefix", 
                      help = "the output prefix")
  parser$add_argument("--method", metavar='', action = "store", 
                      default = "empirical",
    help = paste("either emprical (using the true null) or quantile",
             "(assuming the uniform distribution)"))
  parser$add_argument("--pvalueName", metavar='',
                      action = "store", default = "NA",
    help = paste("the pvalue column name if using the simple quantile method", 
        "assuming the uniform distribution"))
  parser$add_argument("--n", metavar='',
                      action = "store", type = "integer", default = 1000,
    help = "number of replications when sampling from the null distribution")
  parser$add_argument("--setID", metavar='',
                      action = "store", default = "gene",
    help = paste("the column name to uniquely identify each row"))
  parser$add_argument("--outColumns", metavar='',
                      action = "store", default = "chr,gene",
    help = paste("comma separated column names to add to the output,", 
                 "default chr,gene"))
  parser$add_argument("--pattern", metavar='',
                      action = "store", default = "",
    help = paste("specify the column name regular expression patterns forming",
           "the 2x2xm table,", 
          "this is used to select the correct columns to form the 2x2xm table"))
  parser$add_argument("--noLambda",
                      action = "store_true", default = F,
    help = "not to show lambda in the QQ plot")
  parser$add_argument("--nullBoxplot",
                      action = "store_true", default = F,
     help = "show the boxplot of lambdas under the null")
  parser$add_argument("--limit", metavar='',
                      action = "store", type = "double", default = NULL,
                      help = "set the limit of x and y axis of the QQ plot")
  parser$add_argument("--ncore", metavar='',
    action = "store", type = "integer", default = 1,
    help = "number of cores to use")
  parser$add_argument("--alternative", metavar='',
                      action = "store", default = "two.sided",
    help = paste("used in FET or CMH, can be: two.sided, greater, less"))
  parser$add_argument("--excludeIDs", metavar='',
                      action = "store", default = "NA",
        help = paste("a one column file listing the row IDs to be excluded"))
  parser$add_argument("--FDR",
                      action = "store_true", default = F,
                help = "calculate the FDR designed for discrete counts")
  parser$add_argument("--sampling",
                      action = "store", default = "cdf",
      help = "sampling method of null P values: cdf (much faster) or table")

  tryCatch({arguments = parser$parse_args()},
           error = function(e) {
             print(e)
             stop()
           })
  
  inputFile = arguments$inputFile
  outputPrefix = arguments$outputPrefix
  pvalueName = arguments$pvalueName
  method = arguments$method
  nReplication = arguments$n
  columnNamePattern = arguments$pattern
  nolambda = arguments$noLambda
  limit = arguments$limit
  ncore = arguments$ncore
  geneRemove = arguments$excludeIDs
  alternative = arguments$alternative
  nullBoxplot = arguments$nullBoxplot
  returnFDR = arguments$FDR
  sampling = arguments$sampling
  
  setID = arguments$setID
  
  outColumns = arguments$outColumns
  outColumns = strsplit(outColumns, ",")[[1]]
  
  plotFile = paste0(outputPrefix, ".pdf")
  
  input = read.table(inputFile, header = T, as.is = T, comment.char = "")
  
  # make sure setID and outColumns is in the column names
  if (!all(union(setID, outColumns) %in% colnames(input))) {
    print(paste("The column(s)", 
               paste(union(setID, outColumns), collapse = " "), 
                "some column IDs are not in the column names of the input data"))
    stop()
  }
  
  if (geneRemove != "NA") {
    geneRemoveList = as.matrix(read.table(geneRemove, header = F, as.is = T))
    print(geneRemoveList)
    input = input[!(input[, setID] %in% geneRemoveList), ]
  }
  
  pdf(plotFile)
  if (method == "quantile") {
    p = input[, pvalueName]
    p = p[!is.na(p)]
    p[p < 1e-30] = 1e-30
    qqplotQT(p)
  } else if (method == "empirical") {
    columnNames = grep(columnNamePattern, colnames(input), value = T)
    data = as.matrix(input[, columnNames])
    result = qqplotHGAndFDR(data, returnSimulated = F, 
                            nReplication = nReplication,
                            alternative = alternative,
                      nolambda = nolambda, limit = limit, 
                      plotNULL = nullBoxplot,
                      returnFDR = returnFDR,
                      ncore = ncore,
                      sampling = sampling)
    if (returnFDR) {
      FDRResult = result$FDRResult
      out = input[, outColumns, drop = F][result$indexKept, ]
      FDRResult = cbind(out, FDRResult)
      write.table(FDRResult, file = paste0(outputPrefix, ".fdr.tsv"), 
                  row.names = F, col.names = T,
                  sep = "\t", quote = F)
    }
  } else {
    stop(paste("wrong method name:", method))
  }
  dev.off()
}

main()


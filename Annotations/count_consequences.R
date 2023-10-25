tsv <- read.table("output.tsv", header = TRUE, sep = "\t")
tsv <- tsv$REVEL
tsv <- data.frame(tsv)

mayor_65 <- sum(tsv$tsv > 0.65)
cat("Número de variantes con valor mayor de 20:", mayor_65, "\n")

menor_65 <-sum(tsv$tsv >= 0.001 & tsv$tsv <= 0.65)
cat("Número de variantes menor 0.65:", menor_65, "\n")


miss <- read.table("miss.tsv", header = TRUE, sep = "\t")
miss <- miss$CADD_PHRED
miss <- data.frame(miss)

mayor_20 <- sum(miss$miss > 20)
cat("Número de variantes con valor mayor de 20:", mayor_20, "\n")

menor_20 <-sum(miss$miss >= 0.001 & miss$miss <= 20)
cat("Número de variantes menor 20:", menor_20, "\n")

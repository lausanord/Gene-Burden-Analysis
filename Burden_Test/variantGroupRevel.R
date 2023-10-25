pathogenicGroup = function(data, extraParamJason = NULL) {
  # an example of a custom function to define variant sets
  # assuming annovar based annotation, as well as additional annotations 
  # from VEP

  # input
  # data: the data frame where columns are used to define variants of interest
  # extraParamJason: if there are extra parameters or data needed, this 
  #     JASON file provide a way to supply extra parameters.
  
  clinpath <- c("Pathogenic","Pathogenic/Likely_pathogenic","Likely_pathogenic")
  LOF <- c("frameshift_insertion", "frameshift_deletion","stopgain")

  outIndex = logical(dim(data)[1])
  outIndex[ ] = F

  # add clinvar and intervar results
  clinvarIntervalIndex = !grepl("Conflicting",rawData[, "CLNSIG"]) & (
     rawData[, "InterVar_automated"] %in% clinpath | 
     rawData[, "CLNSIG"] %in% clinpath)
  outIndex = outIndex | clinvarIntervalIndex

  outIndex = (data[, "ExonicFunc.refGene"] %in% LOF)

  # add splicing variants
  SpliceAI = max(data$SpliceAI_pred_DS_AG, 
                 data$SpliceAI_pred_DS_AL,
                 data$SpliceAI_pred_DS_DG,
                 data$SpliceAI_pred_DS_DL, na.rm = T)
  dbscSNVScore = max(data$dbscSNV_ADA_SCORE, 
                     data$dbscSNV_RF_SCORE, na.rm = T)
  Splice_HighImpact = (!is.na(SpliceAI) & SpliceAI >= 0.5) | 
                      (!is.na(dbscSNVScore) & dbscSNVScore >= 0.6 |
                        data$LoF == "HC")
  outIndex = outIndex | Splice_HighImpact

  # add missense variants
  outIndex = outIndex | (ExonicFunc.refGene == "nonsynonymous_SNV" & !is.na(REVEL_score) & REVEL_score >= 0.65 - 1e-6)
  
  return(outIndex)
}
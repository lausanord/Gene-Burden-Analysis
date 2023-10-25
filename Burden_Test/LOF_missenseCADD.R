LOF_CADD = function(data, extraParamJason = NULL) {
  # an example of a custom function to define variant sets as LOF, 
  # assuming annovar based annotation
  
  # input
  # data: the data frame where columns are used to define variants of interest
  # extraParamJason: if there are extra parameters or data needed, this 
  #     JASON file provide a way to supply extra parameters.
  
  outIndex = (data[, "ExonicFunc.refGene"] %in% c("nonframeshift_insertion","nonframeshift_deletion",
                                                  "frameshift_deletion", 
                                                  "frameshift_insertion",
                                                  "stopgain")) | 
    (data[, "ExonicFunc.refGene"] == "nonsynonymous_SNV" & (data[, "PHRED_CADD"] == "pathogenic"))
  
  return(outIndex)
}
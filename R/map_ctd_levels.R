map_ctd_levels <- function(results,
                           annot_var="annotLevel"){
  annotLevel <- NULL;
  if(!annot_var %in% names(results)) {
    results[,annotLevel:=ifelse(ctd=="DescartesHuman",2,
                                ifelse(ctd=="HumanCellLandscape",3,NA))]
    # stopper("annot_var not found in results.")
  }
  results[,list(annotLevel=unique(get(annot_var))), keyby="ctd"]
}

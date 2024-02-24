map_ctd_levels <- function(results,
                           annot_var="annotLevel"){
  results[,list(annotLevel=unique(get(annot_var))), keyby="ctd"]
}

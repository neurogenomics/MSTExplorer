order_celltypes <- function(dt,
                            cl=KGExplorer::get_ontology("cl")|>
                              KGExplorer::filter_ontology(
                                keep_descendants = "cell"
                              ),
                            levels=NULL
                            ){
  cl_id <- cl_name <- NULL;
  if(!"cl_id" %in% names(dt)){
    "dt must contain the column 'cl_id'."
  }
  if(!is.null(levels)){
    ## Get celltype order from user-supplied levels
    unmapped_ct <- setdiff(as.character(unique(dt$cl_id)),
                           as.character(levels))
    if(length(unmapped_ct)>0){
      warning(
        paste(
          length(unmapped_ct),"/",length(unique(dt$cl_id)),
          "'cl_id' values in dt could not be mapped",
          "to the user-supplied `levels`."
        )
      )
    }
    levels <- c(levels,unmapped_ct)
  } else {
    ## Get celltype order from dendrogram
    cl_dend <- KGExplorer::ontology_to(cl,
                                       to="dendrogram")
    unmapped_ct <- setdiff(as.character(unique(dt$cl_id)),
                           as.character(labels(cl_dend)))
    if(length(unmapped_ct)>0){
      warning(paste(
        length(unmapped_ct),"/",length(unique(dt$cl_id)),
        "'cl_id' values in dt could not be mapped",
        "to the Cell Ontology."
      ))
    }
    levels <- c(labels(cl_dend),unmapped_ct)
  }
  ## Create ordered factor
  dt[,cl_id:=factor(cl_id,
                    levels=levels,
                    ordered = TRUE)]
  ## order cl_name the same as cl_id
  if("cl_name" %in% names(dt)){
    dt[,cl_name:=factor(
      cl_name,
      levels=unique(dt[order(cl_id),cl_name]),
      ordered = TRUE)]
  }
}

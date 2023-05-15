make_specificity_dt <- function(ctd,
                                annotLevel = 1,
                                shared_genes = NULL,
                                metric = "specificity_quantiles",
                                keep_quantiles = NULL,
                                celltype_col="CellType",
                                verbose = FALSE){
  if(metric=="mean_exp_quantiles" &&
     !"mean_exp_quantiles" %in% names(ctd[[annotLevel]])){
    ctd[[annotLevel]] <- EWCE::bin_specificity_into_quantiles(
      ctdIN =  ctd[[annotLevel]],
      matrix_name = "mean_exp_quantiles",
      numberOfBins = 40,
      verbose = verbose)
  }
  spec <- ctd[[annotLevel]][[metric]]

  value.name <- gsub("s$","",metric)
  spec_dt <- data.table::data.table(as.matrix(spec[shared_genes,]),
                                    keep.rownames = "Gene") |>
    data.table::melt.data.table(id.vars = "Gene",
                                variable.name = celltype_col,
                                value.name = value.name)
  if(!is.null(keep_quantiles)){
    messager("Filtering",metric,v=verbose)
    spec_dt <- spec_dt[get(value.name) %in% keep_quantiles,]
  }
  return(spec_dt)
}

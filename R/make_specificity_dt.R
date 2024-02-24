make_specificity_dt <- function(ctd,
                                annotLevel = 1,
                                shared_genes = NULL,
                                metric = "specificity_quantiles",
                                keep_quantiles = NULL,
                                min_value=NULL,
                                celltype_col="CellType",
                                keep_celltypes = NULL,
                                verbose = FALSE){
  shared_genes <- unique(shared_genes)
  if(metric=="mean_exp_quantiles" &&
     !"mean_exp_quantiles" %in% names(ctd[[annotLevel]])){
    ctd[[annotLevel]] <- EWCE::bin_specificity_into_quantiles(
      ctdIN =  ctd[[annotLevel]],
      matrix_name = "mean_exp_quantiles",
      numberOfBins = 40,
      verbose = verbose)
  }
  spec <- ctd[[annotLevel]][[metric]]
  #### Subset celltypes ####
  keep_celltypes <- keep_celltypes[keep_celltypes %in% colnames(spec)]
  if(length(keep_celltypes)>0){
    messager("Filtering celltypes",v=verbose)
    spec <- spec[,keep_celltypes,drop=FALSE]
  }

  value.name <- gsub("s$","",metric)
  shared_genes <- intersect(shared_genes,
                            rownames(spec))
  spec_dt <- data.table::data.table(as.matrix(spec[shared_genes,,drop=FALSE]),
                                    keep.rownames = "gene_symbol") |>
    data.table::melt.data.table(id.vars = "gene_symbol",
                                variable.name = celltype_col,
                                value.name = value.name)
  if(!is.null(keep_quantiles)){
    messager("Filtering",metric,v=verbose)
    spec_dt <- spec_dt[get(value.name) %in% keep_quantiles,]
  }
  if(!is.null(min_value)){
    messager("Filtering",metric,"by min_value",v=verbose)
    spec_dt <- spec_dt[get(value.name) >= min_value,]
  }
  return(spec_dt)
}

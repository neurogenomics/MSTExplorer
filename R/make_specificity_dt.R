make_specificity_dt <- function(ctd,
                                annotLevel = 1,
                                shared_genes = NULL,
                                metric = "specificity_quantiles"){
  spec <- ctd[[annotLevel]][[metric]]
  data.table::data.table(as.matrix(spec[shared_genes,]),
                                    keep.rownames = "Gene") |>
    data.table::melt.data.table(id.vars = "Gene",
                                variable.name = "celltype_fixed",
                                value.name = gsub("s$","",metric))
}

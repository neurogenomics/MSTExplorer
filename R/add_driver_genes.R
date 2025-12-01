#' Add driver genes
#'
#' This function identifies "driver genes" from phenotype-cell type association analyses.
#' For a given phenotype-cell type pair, driver genes are defined as the intersect of genes that had a
#' phenotype evidence score >0 and were within the top 75th expression specificity percentile
#' (quantiles 30-40 out of 40) for the associated cell type.
#'
#' @param min_value Miniumum specificity quantile to keep.
#' @param metric Which metric to use in the CellTypeDatasets.
#' @param group_var Grouping variable to use when selecting top N genes.
#' Only used when \code{top_n!=NULL}.
#' @param celltype_var The name of the cell type column to merge on.
#' @inheritParams prioritise_targets
#' @inheritParams add_symptom_results
#' @inheritDotParams HPOExplorer::add_genes
#' @export
#' @examples
#' res <- load_example_results()[seq(100)]
#' res <- add_driver_genes(results=res)
add_driver_genes <- function(results = load_example_results(),
                             ctd_list = load_example_ctd(
                               file = paste0("ctd_",unique(results$ctd),".rds"),
                               multi_dataset = TRUE
                               ),
                             annotLevels = map_ctd_levels(results),
                             keep_quantiles = seq(30,40),
                             min_value = NULL,
                             metric = "specificity_quantiles",
                             top_n = NULL,
                             group_var="hpo_id",
                             celltype_var="CellType",
                             ...){
  gene_symbol <- NULL;

  if(metric %in% names(results)){
    messager("specificity_quantiles for driver genes already present in input.")
    return(results)
  }
  results <- HPOExplorer::add_genes(results,
                                    ...)
  #### Find the most cell-type specific genes per cell type per CTD ####
  ctd_dt <- lapply(stats::setNames(names(ctd_list),
                                   names(ctd_list)), function(nm){
                                     make_specificity_dt(
                                       ctd = ctd_list[[nm]],
                                       annotLevel = annotLevels[nm]$annotLevel,
                                       shared_genes = results$gene_symbol,
                                       keep_quantiles = keep_quantiles,
                                       min_value = min_value,
                                       metric = metric)
                                   }) |> data.table::rbindlist(idcol = "ctd")
  #### Find the driver genes for the phenotype-level results ####
  results_ct <- data.table::merge.data.table(x = results,
                                             y = ctd_dt,
                                             by = c("ctd",celltype_var,"gene_symbol"))
  #### Select the top N genes per enrichment result #####
  if(!is.null(top_n)){
    results_ct <- results_ct[,.SD[top_n],
                             by = c(group_var,"ctd",celltype_var)]
  }
  #### Compute number of remaining driver genes per enrichment result ####
  results_ct[,(paste0("n_driver_genes_",group_var)):=data.table::uniqueN(gene_symbol),
             by = c(group_var,"ctd",celltype_var)]
  return(results_ct)
}

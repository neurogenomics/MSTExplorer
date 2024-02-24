add_driver_genes <- function(results = load_example_results(multi_dataset = TRUE),
                             ctd_list = load_example_ctd(
                               file = paste0("ctd_",unique(results$ctd),".rds"),
                               multi_dataset = TRUE),
                             annotLevels = map_ctd_levels(results),
                             keep_quantiles = NULL,
                             min_value = NULL,
                             metric = "specificity_quantiles",
                             top_n = NULL,
                             enrichment_level="hpo_id",
                             ...){
  n_driver_genes <- gene_symbol <- NULL;

  results <- HPOExplorer::add_genes(results,
                                    ...)
  #### Find the most cell-type specific genes per cell type per CTD ####
  ctd_dt <- lapply(stats::setNames(names(ctd_list),
                                   names(ctd_list)), function(nm){
                                     make_specificity_dt(ctd = ctd_list[[nm]],
                                                         annotLevel = annotLevels[nm]$annotLevel,
                                                         shared_genes = results$gene_symbol,
                                                         keep_quantiles = keep_quantiles,
                                                         min_value = min_value,
                                                         metric = metric)
                                   }) |> data.table::rbindlist(idcol = "ctd")
  #### Find the driver genes for the phenotype-level results ####
  results_ct <- data.table::merge.data.table(x = results,
                                             y = ctd_dt,
                                             by = c("ctd","CellType","gene_symbol"))
  #### Select the top N genes per enrichment result #####
  if(!is.null(top_n)){
    results_ct <- results_ct[,.SD[top_n],
                             by = c(enrichment_level,"ctd","CellType")]
  }
  #### Compute number of remaining driver genes per enrichment result ####
  results_ct[,(paste0("n_driver_genes_",enrichment_level)):=data.table::uniqueN(gene_symbol),
             by = c(enrichment_level,"ctd","CellType")]
  return(results_ct)
}

#' Add symptom results
#'
#' Add symptom results to the results data.table.
#' @inheritParams prioritise_targets
#'
#' @export
#' @examples
#' results <- add_symptom_results()
add_symptom_results <- function(results = load_example_results(multi_dataset = TRUE),
                                q_threshold = 0.05,
                                fold_threshold = 2,
                                celltype_col="CellType",
                                ctd_list = load_example_ctd(
                                  file = paste0("ctd_",unique(results$ctd),".rds"),
                                  multi_dataset = TRUE),
                                phenotype_to_genes = HPOExplorer::load_phenotype_to_genes(),
                                annotLevels = list(DescartesHuman=2,
                                                   HumanCellLandscape=3),
                                keep_quantiles = seq(30,40),
                                top_n = NULL,
                                proportion_driver_genes_symptom_threshold=.75
                                ){
  n_genes_hpo_id <- n_genes_disease_id <- n_genes_symptom <- gene_symbol <-
    celltype_symptom <- NULL;

  phenotype_to_genes[,n_genes_hpo_id:=data.table::uniqueN(gene_symbol),
                     by="hpo_id"]
  phenotype_to_genes[,n_genes_disease_id:=data.table::uniqueN(gene_symbol),
                     by="disease_id"]
  phenotype_to_genes[,n_genes_symptom:=data.table::uniqueN(gene_symbol),
                     by=c("hpo_id","disease_id")]
  results <- subset_results(results = results,
                            fold_threshold = fold_threshold,
                            q_threshold = q_threshold)
  if(celltype_col %in% c("cl_id","cl_id")){
    results <- map_celltype(results)
  }
  results <- HPOExplorer::add_genes(results,
                                    phenotype_to_genes = phenotype_to_genes,
                                    allow.cartesian = TRUE)
  results_annot <- data.table::merge.data.table(
    results,
    unique(phenotype_to_genes[,c("hpo_id","disease_id",
                                 "n_genes_hpo_id",
                                 "n_genes_disease_id",
                                 "n_genes_symptom")]),
    by=c("hpo_id","disease_id"))
  #### Add genes that intersect between the
  phenos <- add_driver_genes(results = results_annot,
                             ctd_list = ctd_list,
                             annotLevels = annotLevels,
                             keep_quantiles = keep_quantiles,
                             top_n = top_n,
                             enrichment_level = "hpo_id")
  phenos[,n_driver_genes_symptom:=data.table::uniqueN(gene_symbol),
             by=c("hpo_id","disease_id","ctd","CellType")]
  phenos[,proportion_driver_genes_symptom:=(n_driver_genes_symptom/n_genes_symptom)]
  phenos[proportion_driver_genes_symptom>=proportion_driver_genes_symptom_threshold,
         celltype_symptom:=get(celltype_col)]
  return(phenos)
}

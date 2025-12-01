#' Add symptom results
#'
#' Add symptom results to the results data.table.
#' @param celltype_col Cell type column name in \code{results}.
#' @param annotLevels The annotation level to use within each CTD in
#' \code{ctd_list}.
#' @param keep_quantiles Quantiles to keep in each CellTypeDataset of the
#' \code{ctd_list}.
#' @param proportion_driver_genes_symptom_threshold The minimum proportion of
#' overlap between symptom genes (genes annotated to a phenotype
#' via a specific disease) and the driver genes
#' (genes driving a signficant phenotype-cell type association).
#' @param drop_subthreshold Drop rows that don't meet the
#'  \code{proportion_driver_genes_symptom_threshold} criterion.
#' @inheritParams prioritise_targets
#'
#' @export
#' @examples
#' results <- load_example_results()[seq(5000)]
#' results <- add_symptom_results(results)
add_symptom_results <- function(results = load_example_results(),
                                q_threshold = 0.05,
                                effect_threshold = NULL,
                                celltype_col="CellType",
                                ctd_list = load_example_ctd(
                                  file = paste0("ctd_",unique(results$ctd),".rds"),
                                  multi_dataset = TRUE
                                  ),
                                phenotype_to_genes = HPOExplorer::load_phenotype_to_genes(),
                                annotLevels =  map_ctd_levels(results),
                                keep_quantiles = seq(30,40),
                                top_n = NULL,
                                proportion_driver_genes_symptom_threshold=.25,
                                drop_subthreshold=FALSE
                                ){
  n_genes_hpo_id <- n_genes_disease_id <- n_genes_symptom <- gene_symbol <-
    celltype_symptom <- n_driver_genes_symptom <-
    proportion_driver_genes_symptom <- ..merge_cols <- NULL;

  if("proportion_driver_genes_symptom" %in% names(results)){
    messager("Symptom results already present in input.")
    return(results)
  }
  messager("Adding symptom-level results.")
  phenotype_to_genes[,n_genes_hpo_id:=data.table::uniqueN(gene_symbol),
                     by="hpo_id"]
  phenotype_to_genes[,n_genes_disease_id:=data.table::uniqueN(gene_symbol),
                     by="disease_id"]
  phenotype_to_genes[,n_genes_symptom:=data.table::uniqueN(gene_symbol),
                     by=c("hpo_id","disease_id")]
  results <- subset_results(results = results,
                            effect_threshold = effect_threshold,
                            q_threshold = q_threshold)
  if(celltype_col %in% c("cl_id","cl_name")){
    results <- map_celltype(results)
  }
  results <- HPOExplorer::add_genes(results,
                                    phenotype_to_genes = phenotype_to_genes,
                                    allow.cartesian = TRUE)
  merge_cols <- c("hpo_id",
                  "disease_id",
                  "n_genes_hpo_id",
                  "n_genes_disease_id",
                  "n_genes_symptom")
  if(all(merge_cols %in% colnames(results))){
    results_annot <- results
  }else{
    results_annot <- data.table::merge.data.table(
      results,
      unique(phenotype_to_genes[,..merge_cols]),
      by=c("hpo_id","disease_id"))
  }
  #### Add genes that intersect between the
  phenos <- add_driver_genes(results = results_annot,
                             ctd_list = ctd_list,
                             annotLevels = annotLevels,
                             keep_quantiles = keep_quantiles,
                             top_n = top_n,
                             group_var = "hpo_id")
  phenos[,n_driver_genes_symptom:=data.table::uniqueN(gene_symbol),
             by=c("hpo_id","disease_id","ctd","CellType")]
  phenos[,proportion_driver_genes_symptom:=(n_driver_genes_symptom/n_genes_symptom)]
  phenos[proportion_driver_genes_symptom>=proportion_driver_genes_symptom_threshold,
         celltype_symptom:=get(celltype_col)]
  #### Drop subthreshold symptoms ####
  if(isTRUE(drop_subthreshold) &&
     !is.null(proportion_driver_genes_symptom_threshold)){
    messager("Dropping subthreshold symptoms.")
    phenos <- phenos[
      proportion_driver_genes_symptom>proportion_driver_genes_symptom_threshold,]
  }
  return(phenos)
}

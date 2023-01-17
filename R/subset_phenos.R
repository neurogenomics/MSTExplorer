#' Subset RD EWCE results
#'
#' This subsets  the Rare disease EWCE results by cell type,
#'  q threshold and fold change.
#' @inheritParams ggnetwork_plot_full
#' @inheritParams HPOExplorer::make_phenos_dataframe
#'
#' @returns A data frame of results taken from the main data frame of results.
#' @export
#' @examples
#' phenos <- subset_phenos(cell_type = "Amacrine cells")
subset_phenos <- function(cell_type,
                          results = load_example_results(),
                          hpo = HPOExplorer::get_hpo(),
                          phenotype_to_genes =
                            HPOExplorer::load_phenotype_to_genes(),
                          q_threshold = 0.0005,
                          fold_threshold = 1,
                          verbose = TRUE) {
  # templateR:::source_all()
  # templateR:::args2vars(subset_phenos)
  HPO_ID <- HPO_ID <- HPO_term_valid <- NULL;

  messager("Subsetting results by phenotype.",v=verbose)
  phenos <- get_cell_ontology(cell_type = cell_type,
                              results = results,
                              q_threshold = q_threshold,
                              fold_threshold = fold_threshold,
                              phenotype_to_genes = phenotype_to_genes,
                              hpo = hpo)
  phenos <- phenos[(!is.na(HPO_ID)) & (HPO_term_valid),]
  messager(formatC(nrow(phenos),big.mark = ","),
           "associations remain after filtering.",v=verbose)
  return(phenos)
}

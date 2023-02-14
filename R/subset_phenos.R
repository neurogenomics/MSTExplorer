#' Subset RD EWCE results
#'
#' This subsets  the Rare disease EWCE results by cell type,
#'  q threshold and fold change.
#' @inheritParams ggnetwork_plot_full
#' @inheritParams HPOExplorer::subset_descendants
#' @returns A data frame of results taken from the main data frame of results.
#' @export
#' @importFrom HPOExplorer subset_descendants
#' @examples
#' phenos <- subset_phenos(cell_type = "Amacrine cells",
#'                         ancestor = "Neurodevelopmental delay")
subset_phenos <- function(cell_type,
                          ancestor = NULL,
                          results = load_example_results(),
                          hpo = HPOExplorer::get_hpo(),
                          phenotype_to_genes =
                            HPOExplorer::load_phenotype_to_genes(),
                          q_threshold = 0.0005,
                          fold_threshold = 1,
                          verbose = TRUE) {
  # templateR:::args2vars(subset_phenos)

  #### Subset by q, fold_change, and celltype ####
  phenos <- subset_results(cell_type = cell_type,
                           results = results,
                           q_threshold = q_threshold,
                           fold_threshold = fold_threshold,
                           phenotype_to_genes = phenotype_to_genes,
                           hpo = hpo,
                           verbose = verbose)
  #### Subset by ancestor ####
  phenos <- HPOExplorer::subset_descendants(phenos = phenos,
                                            ancestor = ancestor,
                                            hpo = hpo,
                                            verbose = verbose)
  return(phenos)
}

#' Subset RD EWCE results
#'
#' This subsets  the Rare disease EWCE results by cell type,
#'  q threshold and fold change.
#' @inheritParams ggnetwork_plot_full
#' @inheritParams HPOExplorer::filter_descendants
#' @returns A data frame of results taken from the main data frame of results.
#' @export
#' @importFrom HPOExplorer filter_descendants
#' @examples
#' phenos <- subset_phenos(cell_type = "Cardiomyocytes",
#'                         keep_descendants = "Neurodevelopmental delay")
subset_phenos <- function(cell_type,
                          keep_descendants = NULL,
                          results = load_example_results(),
                          hpo = HPOExplorer::get_hpo(),
                          q_threshold = 0.0005,
                          fold_threshold = 1,
                          verbose = TRUE) {
  #### Subset by q, fold_change, and celltype ####
  phenos <- subset_results(cell_type = cell_type,
                           results = results,
                           q_threshold = q_threshold,
                           fold_threshold = fold_threshold,
                           hpo = hpo,
                           verbose = verbose)
  #### Subset by ancestor ####
  phenos <- HPOExplorer::filter_descendants(phenos = phenos,
                                            keep_descendants = keep_descendants,
                                            hpo = hpo)
  return(phenos)
}

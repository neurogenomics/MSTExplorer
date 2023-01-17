#' Subset EWCE results
#'
#' Subset RD EWCE results data by cell type, fold change and q value.
#' @inheritParams HPOExplorer::make_phenos_dataframe
#' @inheritParams ggnetwork_plot_full
#' @returns A data frame of the selected subset of RD EWCE results
#' with HPO ID column added.
#'
#' @export
#' @importFrom HPOExplorer get_hpo load_phenotype_to_genes add_hpo_id
#' @examples
#' signif_cell_data <- get_cell_ontology(cell_type="Amacrine cells")
get_cell_ontology <- function(cell_type,
                              results = load_example_results(),
                              q_threshold = 0.0005,
                              fold_threshold = 1,
                              hpo = HPOExplorer::get_hpo(),
                              phenotype_to_genes =
                                HPOExplorer::load_phenotype_to_genes()){
  CellType <- fold_change <- NULL;

  message("get_cell_ontology")
  #### Check that celltype is available ####
  if(!any(cell_type %in% unique(results$CellType))){
    cell_orig <- cell_type
    cell_type <- unique(results$CellType)[1]
    message("WARNING!: cell '" ,cell_orig,"' not found in results.\n ",
            "Defaulting to first CellType available: '",cell_type,"'")
  }
  results_sig <- results[(CellType==cell_type) &
                         (q<=q_threshold) &
                         (fold_change>=fold_threshold),]
  #### Check that the table isn't empty after filtering ####
  if(nrow(results_sig)==0){
    stop("ERROR!: phenos table is empty.")
  }
  results_sig <- HPOExplorer::add_hpo_id(
    phenos = results_sig,
    phenotype_to_genes = phenotype_to_genes,
    hpo = hpo)
  return(results_sig)
}

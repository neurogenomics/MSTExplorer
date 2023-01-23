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
#' signif_cell_data <- subset_results(cell_type="Amacrine cells")
subset_results <- function(cell_type,
                           results = load_example_results(),
                           q_threshold = 0.0005,
                           fold_threshold = 1,
                           hpo = HPOExplorer::get_hpo(),
                           phenotype_to_genes =
                             HPOExplorer::load_phenotype_to_genes(),
                           verbose = TRUE){
  CellType <- fold_change <-  HPO_ID <- HPO_ID <- HPO_term_valid <- NULL;

  messager("Subsetting results by q_threshold and fold_change.",v=verbose)
  #### Filter to sig results only ####
  results_sig <- results[(q<=q_threshold) &
                           (fold_change>=fold_threshold),]
  #### Check that celltype is available ####
  if(length(cell_type)==0){
    messager("Skipping cell_type filter.",v=verbose)
  } else if(any(cell_type %in% unique(results_sig$CellType))){
    messager("Subsetting results by cell_type",v=verbose)
    results_sig <- results_sig[CellType %in% cell_type,]
  } else {
    cell_orig <- cell_type
    cell_type <- unique(results_sig$CellType)
    messager("WARNING!: CellType",shQuote(cell_orig),"not found in results.\n ",
            "Defaulting to first CellType available:",shQuote(cell_type))
    results_sig <- results_sig[CellType %in% cell_type,]
  }
  #### Check that the table isn't empty after filtering ####
  if(nrow(results_sig)==0){
    stop("ERROR!: phenos table is empty.")
  }
  #### Add HPO IDs ####
  results_sig <- HPOExplorer::add_hpo_id(
    phenos = results_sig,
    phenotype_to_genes = phenotype_to_genes,
    hpo = hpo,
    verbose = verbose)
  results_sig <- results_sig[(!is.na(HPO_ID)) & (HPO_term_valid),]
  messager(formatC(nrow(results_sig),big.mark = ","),
           "associations remain after filtering.",v=verbose)
  return(results_sig)
}

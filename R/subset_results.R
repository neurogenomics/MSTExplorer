#' Subset EWCE results
#'
#' Subset RD EWCE results data by cell type, fold change and q value.
#' @inheritParams ggnetwork_plot_full
#' @inheritParams KGExplorer::filter_dt
#' @inheritParams HPOExplorer::make_phenos_dataframe
#' @returns A data frame of the selected subset of RD EWCE results
#' with HPO ID column added.
#'
#' @export
#' @examples
#' signif_cell_data <- subset_results(filters=list(CellType = "Cardiomyocytes"))
subset_results <- function(filters = list(cl_name=NULL),
                           results = load_example_results(),
                           q_threshold = 0.0005,
                           effect_threshold = 1,
                           verbose = TRUE){
  effect <-  hpo_id <- hpo_id <- NULL;

  messager("Subsetting results by q_threshold and effect.",v=verbose)
  #### Filter to sig results only ####
  results_sig <- data.table::copy(results)
  if(is.numeric(q_threshold)){
    results_sig <- results_sig[q<q_threshold]
  }
  if(is.numeric(effect_threshold)){
    results_sig <- results_sig[effect>effect_threshold]
  }
  #### Check that celltype is available ####
  if(nrow(results_sig)==0){
    stop("No results remain after filtering")
  } else {
    if(length(unlist(filters))>0){
      messager("Selected",names(filters)[1],":",
               paste("\n -",filters[[1]],collapse=""))
      results_sig2 <- KGExplorer::filter_dt(results_sig,
                                            filters = filters)
      if(nrow(results_sig2)==0){
        cell_orig <- unlist(filters)
        cell_type <- unique(results_sig[[names(filters)[1]]])[[1]]
        messager("WARNING!: CellType",
                 substr(paste(shQuote(cell_orig),collapse = ";"),
                        start = 1, stop = 50),
                 "...",
                 "not found in significant results.\n ",
                 "Defaulting to first CellType available:",shQuote(cell_type))
        filters[[1]] <- cell_type
        results_sig <- KGExplorer::filter_dt(results_sig,
                                             filters = filters)
      } else {
        results_sig <- results_sig2
      }
    }
  }
  #### Check that the table isn't empty after filtering ####
  if(nrow(results_sig)==0){
    stop("ERROR!: phenos table is empty.")
  }
  #### Add HPO IDs ####
  results_sig <- results_sig[(!is.na(hpo_id)) ,]
  messager(formatC(nrow(results_sig),big.mark = ","),
           "associations remain after filtering.",v=verbose)
  return(results_sig)
}

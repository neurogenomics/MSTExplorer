#' merge results
#'
#' merges individual .rds files from the results output directory into one
#' dataframe.
#' @param results_dir The filepath to results .rds files
#' @param description_col The column name for gene list_names (e.g. phenotypes)
#' @return dataframe
#' @export
merge_results <- function(results_dir = "results", description_col = "phenotype") {
  results_merged <- data.frame()
  for (f in list.files(results_dir)) {
    cur = readRDS(paste0(results_dir,f))$results
    descr <-strsplit(f,".rds")
    cur[,description_col] <- descr
    results_merged <- rbind(results_merged,cur)
  }
  return(results_merged)
}

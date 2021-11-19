#' Is not analysed?
#'
#' Takes a gene list name and checks the output results directory to see if that
#' gene list has been analysed yet. The ewce_para function outputs the results
#' for each gene list in to the results directory with name format list_name.rds
#' @param list_name The name of a gene list e.g. "Phenotypic abnormality"
#' @param results_dir The path to results directory
#' @return True or false (bool)
#' @export
is_not_analysed <- function(list_name,results_dir) {
  file_path <- paste0(results_dir,"/",list_name,".rds")
  if (file.exists(file_path)){
    return (FALSE)
  } else {
    return (TRUE)
  }
}

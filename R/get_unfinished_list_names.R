#' Get unfinished list names
#'
#' Gets the unfinished gene list names. It reads the file names in the results
#' directory and uses this to deduce which gene lists have already been analysed.
#' This means you can pause the analysis of multiple gene lists and it will not
#' re-analyse the already completed ones when you start again.
#'
#' @param list_names A char vector of gene list names
#' @param results_dir The directory containing analysed results.
#' @return A character vector of list_names that still need to be analysed.
#' @examples
#' list_names <- HPOExplorer::load_phenotype_to_genes(tempfile())$Phenotype[1:5]
#' results_dir <- tempdir()
#' get_unfinished_list_names(list_names,results_dir)
#'
#' @export
get_unfinished_list_names <- function (list_names, results_dir) {
  list_names_2 = c()
  for (l in list_names) {
    if (is_not_analysed(l,results_dir)) {
      list_names_2 <- append(list_names_2,l)
    }
  }
  return(list_names_2)
}

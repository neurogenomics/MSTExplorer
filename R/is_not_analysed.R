#' Is not analysed?
#'
#' Takes a gene list name and checks the output results directory to see if that
#' gene list has been analysed yet.
#' The \link[MultiEWCE]{ewce_para} function outputs the results
#' for each gene list in to the results directory (\code{save_dir_tmp})
#'  with name format \emph{list_name.rds}.
#' @param list_name The name of a gene list (e.g. "Phenotypic abnormality").
#' @inheritParams ewce_para
#' @returns TRUE or FALSE \code{bool}.

#' @keywords internal
is_not_analysed <- function(list_name,
                            save_dir_tmp) {
  file_path <- make_save_path(save_dir = save_dir_tmp,
                              list_name = list_name)
  return(!file.exists(file_path))
}

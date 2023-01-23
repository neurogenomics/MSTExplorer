#' Is not analysed?
#'
#' Takes a gene list name and checks the output results directory to see if that
#' gene list has been analysed yet.
#' The \link[MultiEWCE]{ewce_para} function outputs the results
#' for each gene list in to the results directory (\code{save_dir_tmp})
#'  with name format \emph{list_name.rds}.
#' @inheritParams ewce_para
#' @returns Vector of TRUE or FALSE.

#' @keywords internal
is_not_analysed <- function(list_names,
                            save_dir_tmp) {
  file_path <- make_save_path(save_dir = save_dir_tmp,
                              list_name = list_names)
  return(!file.exists(file_path))
}

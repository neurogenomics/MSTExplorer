#' Get unfinished list names
#'
#' Gets the unfinished gene list names. It reads the file names in the results
#' directory and uses this to deduce which gene lists
#'  have already been analysed.
#' This means you can pause the analysis of multiple gene lists and it will not
#' re-analyse the already completed ones when you start again.
#'
#' @inheritParams ewce_para
#' @returns A character vector of list_names that still need to be analysed.
#'
#' @export
#' @examples
#' gene_data <- HPOExplorer::load_phenotype_to_genes()
#' list_names <- unique(gene_data$hpo_name)[seq_len(3)]
#' save_dir_tmp <- file.path(tempdir(),"results")
#' ctd <- load_example_ctd()
#' res_files <- ewce_para(ctd = ctd,
#'                        gene_data = gene_data,
#'                        list_names = list_names,
#'                        reps = 10,
#'                        save_dir_tmp = save_dir_tmp)
#' unfinished <- get_unfinished_list_names(list_names = gene_data$hpo_name,
#'                                         save_dir_tmp = save_dir_tmp)
get_unfinished_list_names <- function (list_names,
                                       save_dir_tmp) {
  list_names <- unique(list_names)
  if(is.null(save_dir_tmp)) return(list_names)
  unfinished <- list_names[is_not_analysed(list_names = list_names,
                                            save_dir_tmp = save_dir_tmp)
                          ]
  return(unfinished)
}

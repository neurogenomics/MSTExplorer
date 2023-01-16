#' Get gene list
#'
#' This returns a character vector of genes (a gene list). The gene_data should
#' be a data frame that contains a column of gene list names (e.g. the column
#' may be called "Phenotype"), and a column of genes (e.g. "Gene"). For example,
#' lets call the following data frame \code{phenotype_to_genes}:
#'
#' | Phenotype        | Gene   |
#' | ---------------- | ------ |
#' | "Abnormal heart" | gene X |
#' | "Abnormal heart" | gene Y |
#' | "Poor vision"    | gene Z |
#' | "Poor vision"    | gene Y |
#' @param list_name Names of one or more gene lists.
#' @inheritParams ewce_para
#' @returns A character vector of genes associated with the selected list_name
#' @export
#' @examples
#' gene_data <- HPOExplorer::load_phenotype_to_genes()
#' gene_list <- get_gene_list(list_name="Focal non-motor seizure",
#'                            gene_data = gene_data)
get_gene_list <- function(list_name,
                          gene_data,
                          list_name_column = "Phenotype",
                          gene_column = "Gene"){
  if (!list_name %in% unique(gene_data[[list_name_column]]) ) {
    warning(paste("gene list",
                  list_name, "is not present in" ,
                  list_name_column, "column"))
  }
  return(
    gene_data[[gene_column]][gene_data[[list_name_column]] %in% list_name]
  )
}

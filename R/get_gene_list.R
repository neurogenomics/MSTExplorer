#' get_gene_list
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
#' etc...
#'
#' In the example above, if we wanted to extract genes related to "Abnormal heart",
#' we would run the following example:
#' @examples \dontrun{
#' heart_genes <- get_gene_list(list_name = "Abnormal heart",
#'                              gene_data = phenotype_to_genes,
#'                              list_name_column = "Phenotype",
#'                              gene_column = "Gene")
#' }
#' @param list_name The name of the gene list of interest <string>
#' @param gene_data The data frame of gene list names and associated
#' genes <data.frame>
#' @param list_name_column The name of the column in gene_data that contains
#' gene list names <string>
#' @param gene_column The name of the column in gene_data that contains
#' the genes <string>
#' @returns A charcter vector of genes associated with the selected list_name
#' @export
get_gene_list <- function(list_name,
                          gene_data,
                          list_name_column = "Phenotype",
                          gene_column = "Gene"){
  if (!list_name %in% unique(gene_data[,list_name_column])) {
    warning(paste("gene list", list_name, "is not present in" ,list_name_column, "column"))
  }
  return(paste(gene_data[,gene_column][gene_data[,list_name_column] == list_name]))
}

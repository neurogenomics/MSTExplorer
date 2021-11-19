#' get valid gene lists
#'
#' Gene lists for EWCE must have 4 or more genes. These genes must also be present
#' in the CTD file (see EWCE docs), and must be unique (occasionally a gene list
#' may contain the same gene repeated. This would cause an error).
#'
#' This function gets gene lists that have atleast four unique genes that are also
#' present in the CTD file.
#'
#' @param ctd CTD (cell type data file) see EWCE docs
#' @param list_names Vector of gene list names (character vector)
#' @param gene_data data frame containing a gene column and a column of list_names
#' @param list_name_column The name of the column that contains
#' gene list_names (string)
#' @param gene_column The name of the gene column (e.g. "Gene") (string)
#' @examples \dontrun{
#' remove_invalid_gene_lists(ctd,
#'                           list_names,
#'                           gene_data,
#'                           list_names_column = "Phenotype",
#'                           gene_column = "Gene")
#' }
#' @returns A character vector of list_names that are associated with a valid
#' number of genes
#' @export
get_valid_gene_lists <- function(ctd,
                                 list_names,
                                 gene_data,
                                 list_name_column = "Phenotype",
                                 gene_column = "Gene"){
  ctd_genes = rownames(ctd[[1]]$specificity_quantiles)
  validLists = c()
  for (p in list_names) {
    if (sum(unique(get_gene_list(p,gene_data,list_name_column, gene_column)) %in% ctd_genes) >= 4) {
      validLists = append(validLists, p)
    }
  }
  return(validLists)
}

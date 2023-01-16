#' get valid gene lists
#'
#' Gene lists for EWCE must have 4 or more genes.
#' These genes must also be present
#' in the CTD file (see EWCE docs), and must be unique (occasionally a gene list
#' may contain the same gene repeated. This would cause an error).
#'
#' This function gets gene lists that have at least four unique genes
#' that are also present in the CTD file.
#' @inheritParams ewce_para
#' @inheritParams EWCE::bootstrap_enrichment_test
#' @returns A character vector of list_names that are associated with a valid
#' number of genes
#' @keywords internal
get_valid_gene_lists <- function(ctd,
                                 list_names,
                                 gene_data,
                                 list_name_column = "Phenotype",
                                 gene_column = "Gene",
                                 annotLevel = 1){

  ctd_genes <- rownames(ctd[[annotLevel]]$specificity_quantiles)
  validLists <- lapply(list_names, function(p){
    if (sum(unique(
      get_gene_list(list_name = p,
                    gene_data = gene_data,
                    list_name_column = list_name_column,
                    gene_column = gene_column)
    ) %in% ctd_genes) >= 4) {
      return(p)
    } else {
      return(NULL)
    }
  }) |> unlist()
  return(validLists)
}

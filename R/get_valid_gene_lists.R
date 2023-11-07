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
#'
#' @keywords internal
#' @importFrom data.table uniqueN :=
get_valid_gene_lists <- function(ctd,
                                 list_names,
                                 list_name_column,
                                 gene_data,
                                 annotLevel = 1,
                                 min_genes = 4,
                                 verbose = TRUE){
  gene_symbol <- count <- . <- NULL;

  #### Use secret option to set min_genes in EWCE ####
  Sys.setenv("min_genes"=min_genes)
  #### Check if gene lists are valid ####
  messager("Validating gene lists..",v=verbose)
  ctd_genes <- rownames(ctd[[annotLevel]]$specificity_quantiles)
  shared_genes <- intersect(unique(gene_data$gene_symbol),ctd_genes)
  #### Filter phenotypes #####
  gene_data <- gene_data[get(list_name_column) %in% list_names,
                           ][gene_symbol %in% shared_genes,
                             ]
  gene_counts <- gene_data[,.(count=data.table::uniqueN(gene_symbol)),
                           by=c("hpo_id",list_name_column)]
  validLists <- unique(gene_counts[count>=min_genes,][[list_name_column]])
  gene_data <- gene_data[get(list_name_column) %in% validLists,]
  gene_lists <- lapply(stats::setNames(unique(gene_data[[list_name_column]]),
                                       unique(gene_data[[list_name_column]])),
         function(pheno_i){
      unique(gene_data[get(list_name_column) == pheno_i,]$gene_symbol)
  })
  if(length(gene_lists)==0){
    stopper("No valid gene lists found. ",
            "Check your `gene_data` and `list_names` input.")
  }
  messager(formatC(length(gene_lists),big.mark = ","),"/",
           formatC(length(list_names),big.mark = ","),
           "gene lists are valid.",v=verbose)
  return(gene_lists)
}

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
                                 gene_data,
                                 annotLevel = 1,
                                 verbose = TRUE){
  Gene <- count <- . <- Phenotype <- NULL;

  messager("Validating gene lists..",v=verbose)
  ctd_genes <- rownames(ctd[[annotLevel]]$specificity_quantiles)
  shared_genes <- intersect(unique(gene_data$Gene),ctd_genes)
  gene_data <- gene_data[Phenotype %in% list_names,
                           ][Gene %in% shared_genes,
                             ]
  gene_counts <- gene_data[,.(count=data.table::uniqueN(Gene)),
                           by=c("ID","Phenotype")]
  validLists <- unique(gene_counts[count>=4,]$Phenotype)
  gene_data <- gene_data[Phenotype %in% validLists,]
  gene_lists <- lapply(stats::setNames(unique(gene_data$Phenotype),
                         unique(gene_data$Phenotype)),
         function(pheno_i){
      unique(gene_data[Phenotype == pheno_i,]$Gene)
  })
  messager(formatC(length(gene_lists),big.mark = ","),"/",
           formatC(length(list_names),big.mark = ","),
           "gene lists are valid.",v=verbose)
  return(gene_lists)
}

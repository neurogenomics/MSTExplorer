#' Example prioritised targets
#'
#' @description
#' Example output from the function \link[MultiEWCE]{prioritise_targets}.
#' Used in examples to reduce run time.
#' @source
#' \code{
#' example_targets <- prioritise_targets()
#' usethis::use_data(example_targets, overwrite = TRUE)
#' }
#' @format ontology_index
#' @usage data("example_targets")
"example_targets"


## #' Subphenotype-celltype enrichment results
## #'
## #' @description
## #' Celltype enrichment results for each subphenotype in the HPO annotations,
## #' using the simple \link[MultiEWCE]{gen_overlap} method.
## #' A subphenotype is defined as the combination of an phenotype and
## #'  a particular disease (HPO_ID + LinkID).
## #' @source
## #' \code{
## #' gene_data <- HPOExplorer::load_phenotype_to_genes()
## #' gene_data[,HPO_ID.LinkID:=paste(HPO_ID,LinkID,sep=".")]
## #' subphenotypes <- gen_overlap(gene_data = gene_data,
## #'                              list_name_column = "HPO_ID.LinkID",
## #'                              cores=10,
## #'                              save_dir="~/Downloads")
## #' usethis::use_data(subphenotypes, overwrite = TRUE)
## #' }
## #' @format data.table
## #' @usage data("subphenotypes")
## "subphenotypes"



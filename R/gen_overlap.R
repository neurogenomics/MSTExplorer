#' Generate overlap
#'
#' Compute simply overlap tests for each combination of disease/phenotype
#'  gene setvs. celltype gene set
#'  (determined by the top specificity quantiles for each celltype).
#'
#' \emph{NOTE}:\cr
#' This is a faster but less robust version of \link[MSTExplorer]{gen_results}.
#' It also only requires >=1 gene per disease/phenotype, as opposed to >=4.
#' @param long_format Return results with "union" and "intersection"
#'  genes melted into long format (default: \code{FALSE}).
#'  Otherwise, genes will be collapsed into a list column (\code{TRUE}).
#' @param save_dir Directory to save results to.
#' @inheritParams prioritise_targets
#' @inheritParams ewce_para
#' @returns \link[data.table]{data.table} of all overlap test results.
#'
#' @export
#' @import GeneOverlap
#' @import data.table
#' @importFrom parallel mclapply
#' @examples
#' gene_data <- HPOExplorer::load_phenotype_to_genes()
#' list_names <- unique(gene_data$disease_id)[seq(3)]
#' overlap <- gen_overlap(gene_data = gene_data,
#'                        list_names = list_names)
gen_overlap <- function(gene_data =
                          HPOExplorer::load_phenotype_to_genes(),
                        ctd = load_example_ctd(),
                        list_name_column = "disease_id",
                        gene_column = "gene_symbol",
                        list_names = unique(gene_data[[list_name_column]]),
                        annotLevel = 1,
                        keep_specificity_quantiles=seq(30,40),
                        top_n = NULL,
                        long_format = FALSE,
                        save_dir = tempdir(),
                        force_new = FALSE,
                        cores = 1,
                        verbose = TRUE){
  qval <- pval <- NULL;
  #### Create save path ####
  save_path <- gen_results_save_path(save_dir = save_dir,
                                     prefix = "gen_overlap")
  #### Check if results already exist ####
  if(file.exists(save_path) &&
     isFALSE(force_new)) {
    messager("Results already exist at:",save_path,
             "Use `force_new=TRUE` to overwrite.",v=verbose)
    results_final <- readRDS(save_path)
    return(results_final)
  }
  #### Run new analysis ####
  t1 <- Sys.time()
  ct_genes <- apply(ctd[[annotLevel]]$specificity_quantiles,
                    2,
                    function(x){
                      names(x[x %in% keep_specificity_quantiles])
                    })
  bg <- rownames(ctd[[annotLevel]]$specificity_quantiles)

  #### Iterate overlap tests: data.table ####
  # data.table::setDTthreads(cores)
  # overlap <- gene_data[,gen_overlap_test(ct_genes = ct_genes,
  #                                        list_name = {
  #                                          if(.GRP%%50==0)
  #                                            paste0(.BY," : ",
  #                                                   round(.GRP/.NGRP*100,2),
  #                                                   "%")
  #                                        },
  #                                        dgenes = gene_symbol,
  #                                        long_format = long_format,
  #                                        bg = bg),
  #                      by=list_name_column]

  #### Remove all unnecessary columns to save memory ####
  gene_data <- gene_data[,c(list_name_column,gene_column), with=FALSE]
  #### Subset data to only the list_names ####
  gene_data <- gene_data[get(list_name_column) %in% list_names,]
  messager("Splitting data.",v=verbose)
  split.data.table <- utils::getFromNamespace("split.data.table","data.table")
  gene_data_split <- split.data.table(x = gene_data,
                                      by = list_name_column,
                                      keep.by = FALSE)
  remove(gene_data)
  #### Iterate overlap tests: parallel ####
  fut <- future::plan(future::multisession, workers = cores)
  overlap <- furrr::future_map(.x = gene_data_split,
                     .progress = TRUE,
                     .f = function(.x){
                       gen_overlap_test(ct_genes = ct_genes,
                                        dgenes = unique(.x$gene_symbol),
                                        long_format = long_format,
                                        bg = bg,
                                        verbose = FALSE)
                     })|>
    data.table::rbindlist(use.names = TRUE,
                          idcol = list_name_column)

  #### Multiple testing correction ####
  overlap[,qval:=stats::p.adjust(pval,method = "bonf")]
  #### Get top N ####
  if(!is.null(top_n)){
    overlap <- overlap[,utils::head(.SD, top_n),
                       by = c(list_name_column)]
  }
  #### Report time ####
  messager(difftime(Sys.time(),t1),v = TRUE)
  #### Save results ####
  save_path <- save_results(results = overlap,
                            save_path = save_path,
                            verbose = verbose)
  return(overlap)
}

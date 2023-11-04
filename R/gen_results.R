#' Generate results
#'
#' Generates EWCE results on multiple gene lists in parallel by calling
#' \code{ewce_para}. It allows you to stop the analysis and then continue later
#' from where you left off as it checks the results output directory for
#'finished gene lists and removes them from the input.
#' It also excludes gene lists with
#' less than 4 unique genes (which cause errors in EWCE analysis).
#'
#' The gene_data should be a data frame that contains a column of gene
#' list names (e.g. the column may be called "hpo_name"), and a column of
#' genes (e.g. "gene_symbol"). For example:
#'
#' | hpo_name        | gene_symbol   |
#' | ---------------- | ------ |
#' | "Abnormal heart" | gene X |
#' | "Abnormal heart" | gene Y |
#' | "Poor vision"    | gene Z |
#' | "Poor vision"    | gene Y |
#' etc...
#'
#' For more information on this see docs for get_gene_list
#' (\link[HPOExplorer]{get_gene_lists}).
#' @param save_dir Folder to save merged results in.
#' @param parallel_boot Parallelise at the level of bootstrap iterations,
#' rather than across gene lists.
#' @inheritParams ewce_para
#' @inheritParams EWCE::bootstrap_enrichment_test
#' @returns All results as a dataframe.
#'
#' @export
#' @importFrom HPOExplorer load_phenotype_to_genes
#' @importFrom stringr str_replace_all
#' @examples
#' gene_data <- HPOExplorer::load_phenotype_to_genes()
#' list_names <- unique(gene_data$hpo_id)[seq(5)]
#' ctd <- load_example_ctd()
#' all_results <- gen_results(ctd = ctd,
#'                            gene_data = gene_data,
#'                            list_names = list_names,
#'                            reps = 10)
gen_results <- function(ctd,
                        gene_data,
                        list_name_column = "hpo_id",
                        gene_column = "gene_symbol",
                        list_names = unique(gene_data[[list_name_column]]),
                        bg = unique(gene_data[[gene_column]]),
                        reps = 100,
                        annotLevel = 1,
                        genelistSpecies = "human",
                        sctSpecies = "human",
                        cores = 1,
                        parallel_boot = FALSE,
                        save_dir_tmp = NULL,
                        save_dir = tempdir(),
                        force_new = FALSE,
                        verbose = 1) {

  # devoptera::args2vars(gen_results)

  #### Create save path ####
  save_path <- gen_results_save_path(save_dir = save_dir,
                                     prefix = "gen_results")
  #### Check if results already exist ####
  if(file.exists(save_path) &&
     isFALSE(force_new)) {
    messager("Results already exist at:",save_path,
             "Use `force_new=TRUE` to overwrite.",v=verbose)
    results_final <- readRDS(save_path)
    return(results_final)
  }
  start <- Sys.time()
  #### Run analysis ####
  res_files <- ewce_para(ctd = ctd,
                         list_names = list_names,
                         gene_data = gene_data,
                         list_name_column = list_name_column,
                         gene_column = gene_column,
                         bg = bg,
                         reps = reps,
                         annotLevel= annotLevel,
                         genelistSpecies = genelistSpecies,
                         sctSpecies = sctSpecies,
                         cores = cores,
                         parallel_boot = parallel_boot,
                         save_dir_tmp = save_dir_tmp,
                         force_new = force_new,
                         verbose = verbose)
  #### Merge results into one dataframe ####
  results_final <- merge_results(res_files = res_files,
                                 list_name_column = list_name_column)
  #### Report total time ####
  messager("Done in:",round(difftime(Sys.time(),start,units = "s"), 1),
           "seconds.",v=verbose)
  #### Save merged results ####
  save_path <- save_results(results = results_final,
                            save_path = save_path,
                            verbose = verbose)
  return(results_final)
}

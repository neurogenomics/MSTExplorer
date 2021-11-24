#' Generate results
#'
#' Generates EWCE results on multiple gene lists in parallel by calling
#' \code{ewce_para}. It allows you to stop the analysis and then continue later
#' from where you left off as it checks the results output directory for finished
#' gene lists and removes them from the input. It also excludes gene lists with
#' less than 4 unique genes (which cause errors in ewce analysis).
#'
#' The gene_data should be a data frame that contains a column of gene
#' list names (e.g. the column may be called "Phenotype"), and a column of
#' genes (e.g. "Gene"). For example:
#'
#' | Phenotype        | Gene   |
#' | ---------------- | ------ |
#' | "Abnormal heart" | gene X |
#' | "Abnormal heart" | gene Y |
#' | "Poor vision"    | gene Z |
#' | "Poor vision"    | gene Y |
#' etc...
#'
#' For more information on this see docs for get_gene_list (\code{?get_gene_list})
#'
#' If \code{MergeResutls == TURE}, the function will return a dataframe of all results.
#' No multiple testing corrections are applied to this so it is recommended that
#' they are done after, for example:
#' \code{all_results$q <- stats::p.adjust(all_results$p, method = "BH")}
#'
#' @param ctd The Cell type data file for EWCE analysis (see EWCE docs)
#' @param gene_data The dataframe containing gene list names and associated genes
#' (see docs for get_gene_list for more info).
#' @param list_names The names of each gene list (e.g. "Abnormality of nervous system
#' may be the name of a phenotype assocated gene list)
#' @param background_genes A character vector of background genes (see EWCE docs)
#' @param list_name_column The name of the column in gene_data that contains the
#' gene list names, (e.g. the column may be called "Phenotype" if dealing with
#' phenotype associated gene lists)
#' @param gene_column The name of the column containing genes in the gene_data
#' dataframe. Typically this column is called "Gene"
#' @param results_dir the desired direcory to save results (e.g. "results")
#' @param overwrite_past_analysis overwrite previous results in the results dir \<bool\>
#' @param reps The number of bootstrap reps for EWCE (see ewce docs) \<int\>
#' @param annotLevel The level of cell specificity to select from the CTD,
#' See EWCE docs \<int\>
#' @param genelistSpecies The species ("human"/"mouse") of the gene lists \<string\>
#' @param sctSpecies The species ("human"/"mouse") of the CTD data
#' @param cores The number of cores to run in parallel \<int\>
#' @param MergeResults return merged to single data.frame as a .rds.
#' Note: The function will return merged dataframe even if FALSE \<bool\>
#' @examples
#' gene_data <- HPOExplorer::load_phenotype_to_genes("phenotype_to_genes.txt")
#' ctd <- load_example_CTD()
#' list_names <- unique(gene_data$Phenotype)[1:20]
#' background_genes <- unique(gene_data$Gene)
#' list_name_column <- "Phenotype"
#' gene_column <- "Gene"
#' results_dir <- "results"
#' overwrite_past_analysis <- FALSE
#' MergeResults <- TRUE
#' reps <- 10
#' annotLevel <- 1
#' genelistSpecies <- "human"
#' sctSpecies <- "human"
#' cores <- 1
#'
#' all_results <-MultiEWCE::gen_results(ctd, gene_data, list_names, background_genes,
#'                                      list_name_column, gene_column, results_dir,
#'                                      overwrite_past_analysis, reps, annotLevel,
#'                                      genelistSpecies, sctSpecies, cores,
#'                                      MergeResults)
#'
#'
#' @returns If MergeResults is TRUE, it will return all results as a datframe. If
#' FALSE nothing will be returned, but the individual results will still be saved
#' in the results directory.
#' @export
gen_results <- function(ctd,
                        gene_data,
                        list_names,
                        background_genes,
                        list_name_column = "Phenotype",
                        gene_column = "Gene",
                        results_dir = "results",
                        overwrite_past_analysis = FALSE,
                        reps = 10,
                        annotLevel = 1,
                        genelistSpecies = "human",
                        sctSpecies = "human",
                        cores = 1,
                        MergeResults = FALSE) {

  # remove gene lists that do not have enough valid genes (>= 4)
  list_names <- get_valid_gene_lists(ctd,
                                     list_names = list_names,
                                     gene_data = gene_data,
                                     list_name_column = list_name_column,
                                     gene_column = gene_column)

  # Create results directory and remove finished gene lists
  if (!file.exists(results_dir)) {
    dir.create(results_dir)
  }
  if (!overwrite_past_analysis) {
    list_names <- get_unfinished_list_names(list_names,results_dir)
  }
  # Run analysis
  ewce_para(list_names = list_names,
            gene_data = gene_data,
            list_name_column = "Phenotype",
            gene_column = "Gene",
            results_directory = results_dir,
            ctd_file = ctd,
            background_genes = background_genes,
            bootstrap_reps = reps,
            annotation_Level= annotLevel,
            genes_Species = genelistSpecies,
            ctd_Species = sctSpecies,
            cores = cores)

  # Combine results into a single dataframe
  if (MergeResults) {
    results_final <- merge_results(results_dir = results_dir,
                                   list_name_column = list_name_column)
    saveRDS(results_final,paste0("results_",stringr::str_replace_all(Sys.time(),":","-"),".rds"))
    return(results_final)
  } else {
    results_final <- merge_results(results_dir = results_dir,
                                   list_name_column = list_name_column)
    return(results_final)
  }
}

#' EWCE parallel
#'
#' Runs EWCE in parallel on multiple gene lists.
#' @import parallel
#' @import EWCE
#' @param list_names character vector of gene list names
#' @param gene_lists <list> of character vectors of gene lists, indexed by list_name
#' @param results_directory The desired output filepath for results to be saved
#' @param ctd_file The cell type data object for EWCE analysis (see EWCE docs)
#' @param background_genes The background geneset for EWCE analysis (see EWCE docs)
#' @param bootstrap_reps The number of bootstrap reps <int> (e.g. 100000)
#' @param annotation_Level The level of desired cell resolution from the CTD
#' @param genes_Species The species of gene lists <string> "human" or "mouse"
#' @param ctd_Species "human" or "mouse" <string>
#' @param cores The number of cores to run in parallel (e.g. 8) <int>
#' @return True if analysis was sucessful,
#' saves results at "<results_directory>/<list_name>.rds"
#' @export

ewce_para <- function( list_names,
                       gene_lists,
                       results_directory,
                       ctd_file,
                       background_genes,
                       bootstrap_reps,
                       annotation_Level,
                       genes_Species,
                       ctd_Species,
                       cores) {
  parallel::mclapply(list_names,FUN=function(p,
                                             gene_associations= gene_lists,
                                             results_dir = results_directory,
                                             ctd = ctd_file,
                                             background = background_genes,
                                             reps=bootstrap_reps,
                                             annotLevel = annotation_Level,
                                             genelistSpecies = genes_Species,
                                             sctSpecies = ctd_Species){
    print(p)
    genes = gene_associations[[p]]
    try({
      results = EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                          hits = genes,
                                          bg = background,
                                          reps = reps,
                                          annotLevel = annotLevel,
                                          genelistSpecies=genelistSpecies,
                                          sctSpecies=sctSpecies)
      saveRDS(results, paste0(results_dir,p, ".rds"))
      return(TRUE)})
  },mc.cores=cores)
}

#' EWCE parallel
#'
#' Runs EWCE in parallel on multiple gene lists.
#' @import parallel
#' @import EWCE
#' @param list_names character vector of gene list names
#' @param gene_data data frame of gene list names and genes (see ?get_gene_list)
#' @param list_name_column The name of the gene_data column that has the gene list names
#' @param gene_column The name of the gene_data column that contains the genes
#' @param results_directory The desired output filepath for results to be saved
#' @param ctd_file The cell type data object for EWCE analysis (see EWCE docs)
#' @param background_genes The background geneset for EWCE analysis (see EWCE docs)
#' @param bootstrap_reps The number of bootstrap reps `int` (e.g. 100000)
#' @param annotation_Level The level of desired cell resolution from the CTD
#' @param genes_Species The species of gene lists `string` "human" or "mouse"
#' @param ctd_Species "human" or "mouse" `string`
#' @param cores The number of cores to run in parallel (e.g. 8) `int`
#' @return True if analysis was sucessful. The results will then be saved
#' at "(results_directory)/(list_name).rds"
#'
#' @examples
#' gene_data <- HPOExplorer::load_phenotype_to_genes("phenotype_to_genes.txt")
#' ctd <- MultiEWCE::load_example_CTD()
#' list_names <- unique(gene_data$Phenotype)[1:10]
#' results <- ewce_para(list_names, gene_data, results_directory ="results",ctd_file = ctd,
#'           background_genes = unique(gene_data$Gene), bootstrap_reps = 10,
#'           annotation_Level = 1, genes_Species = "human", ctd_Species = "human",
#'           cores = 1)
#'
#'
#' @export
ewce_para <- function( list_names,
                       gene_data,
                       list_name_column = "Phenotype",
                       gene_column = "Gene",
                       results_directory,
                       ctd_file,
                       background_genes,
                       bootstrap_reps,
                       annotation_Level,
                       genes_Species,
                       ctd_Species,
                       cores) {
  parallel::mclapply(list_names,FUN=function(p,
                                             gene_associations= gene_data,
                                             lst_nm_col = list_name_column,
                                             gn_col = gene_column,
                                             results_dir = results_directory,
                                             ctd = ctd_file,
                                             background = background_genes,
                                             reps=bootstrap_reps,
                                             annotLevel = annotation_Level,
                                             genelistSpecies = genes_Species,
                                             sctSpecies = ctd_Species){
    print(p)
    genes = get_gene_list(p,gene_associations,lst_nm_col, gn_col)
    try({
      results = EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                          hits = genes,
                                          bg = background,
                                          reps = reps,
                                          annotLevel = annotLevel,
                                          genelistSpecies=genelistSpecies,
                                          sctSpecies=sctSpecies)
      saveRDS(results, paste0(results_dir,"/",p, ".rds"))
      return(TRUE)})
  },mc.cores=cores)
}

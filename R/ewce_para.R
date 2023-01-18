#' EWCE parallel
#'
#' Runs EWCE in parallel on multiple gene lists.
#' @param ctd Cell Type Data List generated using
#'  \link[EWCE]{generate_celltype_data}.
#' @param list_names character vector of gene list names.
#' @param gene_data data frame of gene list names and genes
#' (see \link[HPOExplorer]{get_gene_lists}).
#' @param list_name_column The name of the gene_data column
#' that has the gene list names.
#' @param gene_column The name of the gene_data column that contains the genes.
#' @param save_dir_tmp Folder to save intermediate results files to
#' (one file per gene list). Set to \code{NULL} to skip saving temporary files.
#' @param cores The number of cores to run in parallel (e.g. 8) \code{int}.
#' @inheritParams EWCE::bootstrap_enrichment_test
#' @returns Paths to saved results at "(save_dir)/(list_name).rds"
#' (when \code{save_dir!=NULL}), or a nested list of results
#' (when \code{save_dir==NULL}).
#'
#' @export
#' @importFrom HPOExplorer load_phenotype_to_genes get_gene_lists
#' @importFrom stats setNames
#' @importFrom parallel mclapply
#' @importFrom EWCE bootstrap_enrichment_test
#' @examples
#' gene_data <- HPOExplorer::load_phenotype_to_genes()
#' ctd <- MultiEWCE::load_example_ctd()
#' list_names <- unique(gene_data$Phenotype)[seq_len(3)]
#' res_files <- ewce_para(ctd = ctd,
#'                        gene_data = gene_data,
#'                        list_names = list_names,
#'                        reps = 10)
ewce_para <- function(ctd,
                      gene_data,
                      list_name_column = "Phenotype",
                      gene_column = "Gene",
                      list_names = unique(gene_data[[list_name_column]]),
                      bg = unique(gene_data[[gene_column]]),
                      reps=100,
                      annotLevel=1,
                      genelistSpecies="human",
                      sctSpecies="human",
                      save_dir_tmp = tempdir(),
                      cores=1,
                      verbose=FALSE) {
  # templateR:::source_all()
  # templateR:::args2vars(ewce_para)

  if(!is.null(save_dir_tmp)){
    dir.create(save_dir_tmp, showWarnings = FALSE, recursive = TRUE)
  }
  list_names <- unique(list_names)
  res_files <- parallel::mclapply(stats::setNames(list_names,
                                                  list_names),
                     FUN=function(p){
    i <- which(list_names==p)
    message_parallel("Analysing: ",p," (",i,"/",length(list_names),")")
    genes <- HPOExplorer::get_gene_lists(
      phenotypes = p,
      phenotype_to_genes = gene_data)[["Gene"]]
    tryCatch({
      results <- EWCE::bootstrap_enrichment_test(
        sct_data = ctd,
        hits = genes,
        bg = bg,
        reps = reps,
        annotLevel = annotLevel,
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies,
        verbose = verbose>1 )
      if(!is.null(save_dir_tmp)){
        save_path <- make_save_path(save_dir = save_dir_tmp,
                                    list_name = p)
        saveRDS(results, save_path)
        return(save_path)
      } else {
        return(results)
      }
      },
      error=function(e){message(e);NULL})
  },mc.cores=cores)
  return(res_files)
}

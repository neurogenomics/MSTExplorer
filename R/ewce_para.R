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
#' @param force_new Overwrite previous results
#'  in the \code{save_dir_tmp}.
#' @param cores The number of cores to run in parallel (e.g. 8) \code{int}.
#' @param min_genes Minimum number of genes per list (default: 4)
#' @inheritParams gen_results
#' @inheritParams EWCE::bootstrap_enrichment_test
#' @inheritDotParams EWCE::bootstrap_enrichment_test
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
#' list_names <- unique(gene_data$hpo_id)[seq(3)]
#' res_files <- ewce_para(ctd = ctd,
#'                        gene_data = gene_data,
#'                        list_names = list_names,
#'                        reps = 10)
ewce_para <- function(ctd,
                      gene_data,
                      list_name_column = "hpo_id",
                      gene_column = "gene_symbol",
                      list_names = unique(gene_data[[list_name_column]]),
                      reps=100,
                      annotLevel=1,
                      force_new=FALSE,
                      genelistSpecies="human",
                      sctSpecies="human",
                      bg = get_bg(species1 = genelistSpecies,
                                  species2 = sctSpecies,
                                  overwrite = force_new),
                      min_genes = 4,
                      save_dir_tmp = tempdir(),
                      parallel_boot = FALSE,
                      cores = 1,
                      verbose = FALSE,
                      ...) {
  # devoptera::args2vars(ewce_para)

  if(!is.null(save_dir_tmp)){
    dir.create(save_dir_tmp, showWarnings = FALSE, recursive = TRUE)
  }
  #### remove gene lists that do not have enough valid genes (>= 4) ####
  gene_lists <- get_valid_gene_lists(ctd = ctd,
                                     annotLevel = annotLevel,
                                     list_names =  unique(list_names),
                                     list_name_column = list_name_column,
                                     gene_data = gene_data,
                                     min_genes = min_genes,
                                     verbose = verbose)

  #### Create results directory and remove finished gene lists ####
  list_names <- if (isFALSE(force_new)) {
    get_unfinished_list_names(list_names = names(gene_lists),
                              save_dir_tmp = save_dir_tmp)
  } else {
    names(gene_lists)
  }
  #### Parallelise at different levels ####
  if(isTRUE(parallel_boot)){
    no_cores <- cores
    cores <- 1
  } else {
    no_cores <- 1
  }
  #### Report on background
  messager("Background contains",formatC(length(bg),big.mark = ","),"genes.",
           v=verbose)
  #### Iterate EWCE ####
  res_files <- parallel::mclapply(stats::setNames(list_names,
                                                  list_names),
                     FUN=function(p){
    i <- which(list_names==p)
    genes <- gene_lists[[p]]
    messager("Analysing: ",shQuote(p),
                     " (",i,"/",length(list_names),"): ",
                     formatC(length(genes),big.mark = ",")," genes.",
             parallel = TRUE)
    tryCatch({
      results <- EWCE::bootstrap_enrichment_test(
        sct_data = ctd,
        hits = genes,
        bg = bg,
        reps = reps,
        annotLevel = annotLevel,
        genelistSpecies = genelistSpecies,
        sctSpecies = sctSpecies,
        no_cores = no_cores,
        verbose = verbose>1,
        ...)
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

#' Generate overlap
#'
#' Compute simply overlap tests for each combination of disease/phenotype
#'  gene setvs. celltype gene set
#'  (determined by the top specificity quantiles for each celltype).
#'
#' \emph{NOTE}:\cr
#' This is a faster but less robust version of \link[MultiEWCE]{gen_results}.
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
#' list_names <- unique(gene_data$LinkID)[seq_len(3)]
#' overlap <- gen_overlap(gene_data = gene_data,
#'                        list_names = list_names)
gen_overlap <- function(gene_data =
                          HPOExplorer::load_phenotype_to_genes(),
                        ctd = load_example_ctd(),
                        list_name_column = "LinkID",
                        gene_column = "Gene",
                        list_names = unique(gene_data[[list_name_column]]),
                        annotLevel = 1,
                        keep_specificity_quantiles=seq(30,40),
                        top_n = NULL,
                        long_format = FALSE,
                        save_dir = tempdir(),
                        cores = 1,
                        verbose = TRUE){
  LinkID <- qval <- pval <-  NULL;

  ct_genes <- apply(ctd[[annotLevel]]$specificity_quantiles,
                    2,
                    function(x){
                      x[x %in% keep_specificity_quantiles]
                    })
  bg <- rownames(ctd[[annotLevel]]$specificity_quantiles)
  #### Iterate overlap tests ####
  overlap <- parallel::mclapply(stats::setNames(list_names,
                                                list_names),
                                function(d){
    messager("Testing overlap: ",d, parallel = TRUE)
    d2 <- gene_data[LinkID==d]
    dgenes <- unique(d2$Gene)

    lapply(ct_genes, function(cgenes){
      r <- GeneOverlap::newGeneOverlap(
        listA = names(cgenes),
        listB = dgenes,
        genome.size = length(bg)) |>
        GeneOverlap::testGeneOverlap()
      data.table::data.table(
        intersection=if(long_format) r@intersection else list(r@intersection),
        intersection_size=length(r@intersection),
        union=if(long_format) r@union else list(r@union),
        union_size=length(r@union),
        pval=r@pval,
        odds.ratio=r@odds.ratio,
        jaccard=r@Jaccard)
    }) |>
      data.table::rbindlist(use.names = TRUE, idcol = "CellType")
  }, mc.cores = cores) |>
    data.table::rbindlist(use.names = TRUE, idcol = "LinkID")
  #### Multiple testing correction ####
  overlap[,qval:=stats::p.adjust(pval,method = "bonf")]
  #### Get top N ####
  if(!is.null(top_n)){
    overlap <- overlap[,utils::head(.SD, top_n),
                       by = c(list_name_column)]
  }
  save_path <- save_results(results = overlap,
                            save_dir = save_dir,
                            prefix = "gen_overlap_",
                            verbose = verbose)
  return(overlap)
}

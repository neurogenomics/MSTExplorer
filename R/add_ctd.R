#' Add CellTypeDataset
#'
#' Annotate genes in \code{results} with data from a CellTypeDataset (CTD).
#' The following columns will be added:
#' \itemize{
#' \item{"specificity"} Cell-type specificity score.
#' \item{"specificity_quantiles"} Cell-type specificity quantile.
#' \item{"mean_exp"} Mean expression per cell-type.
#' \item{"mean_exp_quantiles"} Mean expression quantile per cell-type.
#' }
#' @param keep_celltypes A character vector of cell types to keep.
#' @inheritParams prioritise_targets
#' @inheritParams plot_report
#' @inheritParams data.table::merge.data.table
#'
#' @export
#' @examples
#' results <- load_example_results()[seq(100),]
#' results2 <- add_ctd(results=results)
add_ctd <- function(results=load_example_results(),
                    ctd=load_example_ctd(),
                    annotLevel=length(ctd),
                    keep_specificity_quantiles=NULL,
                    keep_mean_exp_quantiles=NULL,
                    keep_celltypes=NULL,
                    rep_dt=NULL,
                    all.x=FALSE,
                    allow.cartesian = TRUE,
                    by=c("gene_symbol","CellType"),
                    verbose=TRUE){
  gene_symbol <- NULL;

  results <- HPOExplorer::add_genes(phenos = results,
                                    all.x = all.x,
                                    allow.cartesian = allow.cartesian)
  #### Identify genes within CTD ####
  shared_genes <- intersect(results$gene_symbol,
                            rownames(ctd[[annotLevel]]$specificity))
  #### Merge genes with phenotype/celltype results ####
  results <- results[gene_symbol %in% shared_genes,]
  #### Format CTD data ####
  spec_df <- make_specificity_dt(ctd = ctd,
                                 annotLevel = annotLevel,
                                 shared_genes = shared_genes,
                                 keep_celltypes = keep_celltypes,
                                 metric = "specificity")
  specq_df <- make_specificity_dt(ctd = ctd,
                                  annotLevel = annotLevel,
                                  shared_genes = shared_genes,
                                  keep_quantiles = keep_specificity_quantiles,
                                  keep_celltypes = keep_celltypes,
                                  metric = "specificity_quantiles")
  exp_df <- make_specificity_dt(ctd = ctd,
                                annotLevel = annotLevel,
                                shared_genes = shared_genes,
                                keep_celltypes = keep_celltypes,
                                metric = "mean_exp")
  expq_df <- make_specificity_dt(ctd = ctd,
                                 annotLevel = annotLevel,
                                 shared_genes = shared_genes,
                                 keep_quantiles = keep_mean_exp_quantiles,
                                 keep_celltypes = keep_celltypes,
                                 metric = "mean_exp_quantiles")
  #### Merge: specificity ####
  results <- results |>
    data.table::merge.data.table(y = spec_df,
                                 by = by,
                                 all.x = all.x,
                                 allow.cartesian = allow.cartesian)
  #### Merge: specificity_quantiles ####
  results <- results |>
    data.table::merge.data.table(y = specq_df,
                                 by = by,
                                 all.x = all.x,
                                 allow.cartesian = allow.cartesian)
  if(!is.null(rep_dt)){
    rep_dt <- report(dt = results,
                     rep_dt = rep_dt,
                     step = "keep_specificity_quantiles",
                     verbose = verbose)
  }
  #### Merge: mean_exp ####
  results <- results |>
    data.table::merge.data.table(y = exp_df,
                                 by = by)
  #### Merge: mean_exp_quantiles ####
  results <- results |>
    data.table::merge.data.table(y = expq_df,
                                 by = by)
  if(!is.null(rep_dt)){
    rep_dt <- report(dt = results,
                     rep_dt = rep_dt,
                     step = "keep_mean_exp_quantiles",
                     verbose = verbose)
  }
  #### Return ####
  return(
    list(results=results,
         rep_dt=rep_dt,
         shared_genes=shared_genes)
  )
}

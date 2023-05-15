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
#' @inheritParams prioritise_targets
#' @inheritParams data.table::merge.data.table
#'
#' @export
#' @examples
#' res <- load_example_results("Descartes_All_Results_extras.rds")
#' res2 <- add_ctd(results=res)
add_ctd <- function(results=load_example_results(),
                    ctd=load_example_ctd(),
                    annotLevel=length(ctd),
                    keep_specificity_quantiles=NULL,
                    keep_mean_exp_quantiles=NULL,
                    rep_dt=NULL,
                    all.x=FALSE,
                    by=c("Gene","CellType"),
                    verbose=TRUE){
  # devoptera::args2vars(add_ctd)
  Gene <- NULL;

  #### Identify genes within CTD ####
  shared_genes <- intersect(results$Gene,
                            rownames(ctd[[annotLevel]]$specificity))
  #### Merge genes with phenotype/celltype results ####
  results <- results[Gene %in% shared_genes,]
  #### Format CTD data ####
  spec_df <- make_specificity_dt(ctd = ctd,
                                 annotLevel = annotLevel,
                                 shared_genes = shared_genes,
                                 metric = "specificity")
  specq_df <- make_specificity_dt(ctd = ctd,
                                  annotLevel = annotLevel,
                                  shared_genes = shared_genes,
                                  keep_quantiles = keep_specificity_quantiles,
                                  metric = "specificity_quantiles")
  exp_df <- make_specificity_dt(ctd = ctd,
                                annotLevel = annotLevel,
                                shared_genes = shared_genes,
                                metric = "mean_exp")
  expq_df <- make_specificity_dt(ctd = ctd,
                                 annotLevel = annotLevel,
                                 shared_genes = shared_genes,
                                 keep_quantiles = keep_mean_exp_quantiles,
                                 metric = "mean_exp_quantiles")
  #### Merge: specificity ####
  results <- results |>
    data.table::merge.data.table(y = spec_df,
                                 by = by,
                                 all.x = all.x)
  #### Merge: specificity_quantiles ####
  results <- results |>
    data.table::merge.data.table(y = specq_df,
                                 by = by,
                                 all.x = all.x)
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

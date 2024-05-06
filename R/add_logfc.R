#' Add log fold-change
#'
#' Compute \code{log2(fold-change)} based on the effect size.
#' @param force_new Force the creation of a new \code{logFC} column
#'  even when one already exists.
#' @inheritParams prioritise_targets
#' @export
#' @examples
#' results = load_example_results()[seq(5000)]
#' add_logfc(results)
add_logfc <- function(results,
                      effect_var="estimate",
                      force_new=FALSE){
  logFC <- NULL;
  if(!"logFC" %in% names(results) || isTRUE(force_new)){
    messager("Adding logFC column.")
    # results[,logFC:=log10(scales::rescale(estimate,c(.Machine$double.xmin, 1)))]
    # results[,logFC:=logFC/mean(results$logFC)]
    results[,logFC:=(get(effect_var)+abs(min(get(effect_var))))]
    results[,logFC:=(log2(logFC/mean(results$logFC)))]
    # hist(results$logFC)
  } else {
    messager("logFC already exists in results.",
             "Use `force_new=TRUE` to overwrite.")
  }
}

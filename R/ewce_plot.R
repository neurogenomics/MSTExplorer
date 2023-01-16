#' EWCE plot
#'
#' A shallow wrapper for the \link[EWCE]{ewce_plot} function.
#' @inheritParams EWCE::ewce_plot
#' @inheritDotParams EWCE::ewce_plot
#' @returns A bar chart with dendrogram of EWCE results in each cell type.
#'
#' @export
#' @importFrom EWCE ewce_plot
#' @examples
#' full_results <- EWCE::example_bootstrap_results()
#' plt <- MultiEWCE::ewce_plot(total_res = full_results$results)
ewce_plot <-function (total_res,
                      mtc_method = "bonferroni",
                      ctd = NULL,
                      ...) {
  requireNamespace("ggplot2")

  EWCE::ewce_plot(total_res = total_res,
                  mtc_method = mtc_method,
                  ctd = ctd,
                  ...)
}

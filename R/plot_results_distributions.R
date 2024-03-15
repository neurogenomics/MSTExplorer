#' Plot results distributions
#'
#' Plot the distributions of the test statistics in the phenotype-cell type
#' association results.
#' @param log10_transform log10 transform \code{x_var}.
#' @inheritParams plot_
#' @inheritParams prioritise_targets
#' @inheritParams ggplot2::facet_wrap
#' @inheritParams ggplot2::geom_histogram
#' @inheritParams KGExplorer::plot_save
#' @inheritDotParams KGExplorer::plot_save
#' @export
#' @examples
#' results <- load_example_results()
#' ct <- unique(results$CellType)[seq(9)]
#' results <- results[CellType %in% ct]
#' results[,estimate_scaled:=scales::rescale(estimate,c(1,5))]
#' out <- plot_results_distributions(results, x_var="estimate_scaled")
plot_results_distributions <- function(results = load_example_results(),
                                       x_var = "estimate",
                                       facets = "CellType",
                                       binwidth = NULL,
                                       scales = "fixed",
                                       log10_transform=TRUE,
                                       show_plot = TRUE,
                                       save_path=NULL,
                                       ...){

 plt <- ggplot2::ggplot(results,
                 ggplot2::aes(x=!!ggplot2::sym(x_var)
                              )
                 ) +
    ggplot2::geom_histogram(binwidth = binwidth) +
    ggplot2::facet_wrap(facets = facets,
                        scales = scales) +
    ggplot2::theme_bw()
  if(log10_transform){
    plt <- plt + ggplot2::scale_x_log10()
  }
 if(show_plot) methods::show(plt)
 KGExplorer::plot_save(plt = plt,
                       save_path = save_path,
                       ...)
 return(plt)
}

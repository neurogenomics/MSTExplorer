plot_density_cor <- function(res,
                             x="p_HumanCellLandscape",
                             y="p_DescartesHuman",
                             point.args = list(alpha=.1),
                             density_alpha=.7,
                             min_colors=15,
                             log_vars=FALSE){
  res <- res[!is.na(get(x)) & !is.na(get(y))]
  p<- res |>
    ggstatsplot::ggscatterstats(x=!!ggplot2::sym(x),
                                y=!!ggplot2::sym(y),
                                xsidehistogram.args = list(fill="magenta", color="white"),
                                ysidehistogram.args = list(fill="blue", color="white"),
                                point.args = point.args) +
    ggplot2::geom_density_2d_filled(alpha=density_alpha) +
    ggplot2::scale_fill_manual(
      values = c(ggplot2::alpha("white",0),
                 ggplot2::alpha(pals::gnuplot(min_colors)[-1],density_alpha)
    ))
  if(isTRUE(log_vars)){
    p <- p+
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10()
  }
  return(p)
}

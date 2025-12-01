plot_density_cor <- function(res,
                             x="p_HumanCellLandscape",
                             y="p_DescartesHuman",
                             agg_var=NULL,
                             add_density=TRUE,
                             rasterize_points=TRUE,
                             point.args = list(alpha=.5),
                             dpi = 100,
                             density_alpha=.7,
                             marginal=FALSE,
                             min_colors=15,
                             log_vars=FALSE,
                             downsample=NULL,
                             stats_idx=c(1,3,4,7),
                             show.legend=FALSE
                             ){
  . <- NULL;

  res <- data.table::copy(res)
  res <- res[!is.na(get(x)) & !is.na(get(y))]


  x_orig <- x
  y_orig <- y

  ### Downsample for plotting
  if (!is.null(downsample) && downsample < nrow(res)) {
    downsample <- min(downsample, nrow(res))
    messager("Downsampling to", downsample, "points.")
    res <- res[sample(.N, downsample)]
  }
  res <- res[is.finite(get(x)) & is.finite(get(y))]


  # Aggregate to reduce plot size
  if(!is.null(agg_var)){
    agg_dt <- res[,.(x=mean(get(x),na.rm=TRUE),
                     y=mean(get(y),na.rm=TRUE)),
                  by=agg_var]
    res <- agg_dt
    x <- "x"
    y <- "y"
  }

  p <- res |>
    ggstatsplot::ggscatterstats(x=!!ggplot2::sym(x),
                                y=!!ggplot2::sym(y),
                                marginal=marginal,
                                marginal.type = "histogram",
                                xsidehistogram.args = list(fill="darkgrey", color="white"),
                                ysidehistogram.args = list(fill="darkgrey", color="white"),
                                bf.message=FALSE,
                                point.args = point.args
                                )  +
    ggplot2::labs(title=paste(x_orig,"vs.",y_orig),
                  x=x_orig, y=y_orig,
                  fill=if(add_density) "Density")

  if(isTRUE(add_density)){
    p <- p + ggplot2::geom_density_2d_filled(alpha=density_alpha,
                                    show.legend = show.legend) +
      ggplot2::scale_fill_manual(
        values = c(ggplot2::alpha("white",0),
                   ggplot2::alpha(pals::gnuplot(min_colors)[-1],density_alpha)
        ))
  }

  p <- filter_ggstatsplot_subtitle(p, stats_idx = stats_idx)

  if (rasterize_points){
    p <- ggrastr::rasterize(p, layers = "Point", dpi = dpi)
  }


  if(isTRUE(log_vars)){
    p <- p+
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10()
  }
  return(p)
}

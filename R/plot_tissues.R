plot_tissues <- function(results,
                         fill_var="n_enrichments",
                         facet_var=NULL,
                         by=c("cl_name",facet_var),
                         label_size=2,
                         y_adjust=1,
                         hjust = -0.2,
                         types=c("tile","bar")){
  cl_name <- uberon_ancestor_name <- n_enrichments <- NULL;
  res <- list()
  {
    if(("tile" %in% types) ||
        any(c(fill_var,facet_var,by) %in% c("uberon_name","uberon_id"))){
      results <- map_tissue(results = results)
    }
    if(fill_var!="n_enrichments"){
      by <- c(by,fill_var)
    }
    by <- intersect(by,names(results))
    plot_dat <- results[q<0.05, list(n_enrichments=.N), by=by]
    res[["data"]] <- plot_dat
  }

  if("tile" %in% types){
    res[["tile_plot"]] <- ggplot2::ggplot(plot_dat,
                              ggplot2::aes(x=cl_name,
                                           y=uberon_ancestor_name,
                                           fill=!!ggplot2::sym(fill_var)
                              )
    ) +
      ggplot2::geom_tile(show.legend = FALSE) +
      ggplot2::scale_fill_viridis_c(option = "plasma",end = .95) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 1, vjust=0.5),
                     strip.background = ggplot2::element_rect(fill="transparent")
      )
    #### Add facets ####
    if(!is.null(facet_var)){
      res[["tile_plot"]]  <- res[["tile_plot"]] +
        ggplot2::facet_wrap(facets=facet_var)
    }
  }

  if("bar" %in% types){
    res[["bar_plot"]] <- ggplot2::ggplot(plot_dat,
                             ggplot2::aes(x=cl_name,
                                          y=n_enrichments,
                                          fill=!!ggplot2::sym(fill_var)
                             )
    ) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::geom_text(ggplot2::aes(label = n_enrichments,
                                      y = n_enrichments*y_adjust,
                                      group=cl_name,
                                      fill=NULL),
                         hjust = -0.2,
                         # label.padding = ggplot2::unit(0.1,"lines"),
                         angle = 90,
                         size = label_size) +
      ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult = c(0,.2))) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 1,
                                                         vjust=0.5),
                     axis.ticks.x = ggplot2::element_blank(),
                     strip.background = ggplot2::element_rect(fill = "transparent")
      )
    #### Add facets ####
    if(!is.null(facet_var)){
      res[["bar_plot"]] <- res[["bar_plot"]] +
        ggplot2::facet_wrap(facets=facet_var)
    }
  }

  #### Return ####
  return(res)
}

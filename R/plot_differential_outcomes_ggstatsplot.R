plot_differential_outcomes_ggstatsplot <- function(plot_dat,
                                                   remove_facets,
                                                   x_var,
                                                   y_var,
                                                   facet_var,
                                                   max_facets,

                                                   type = "nonparametric",
                                                   centrality.type = type,
                                                   centrality.plotting = FALSE,
                                                   pairwise.comparisons = TRUE,
                                                   ...){
  requireNamespace("paletteer")
  requireNamespace("ggstatsplot")
  median_y <- pct <- N <- NULL;
  # valid_facets <- plot_dat[, list(N=.N,
  #                                 sd=sd(get(y_var))
  #                                 ),
  #                          by=c(facet_var,x_var)
  #                      ][,list(valid=all(N>1& sd>0)),
  #                              by=facet_var
  #                        ][valid==TRUE][[facet_var]]
  plot_dat <- plot_dat[!is.na(get(x_var)) & !is.na(get(y_var)) & !is.na(get(facet_var))]
  valid_facets <- unique(plot_dat[[facet_var]])
  messager("Valid facets:",length(valid_facets))
  valid_dat <- limit_facets(dat = plot_dat,
                            facet_var = facet_var,
                            facet_subset = valid_facets,
                            max_facets = max_facets)
  valid_pal <- subset(paletteer::palettes_d_names,
                      length>=length(unique(plot_dat[[x_var]])))|>
    data.frame()
  valid_dat <- valid_dat[!get(facet_var) %in% remove_facets]
  plts_list <- lapply(stats::setNames(unique(valid_dat[[facet_var]]),
                                      unique(valid_dat[[facet_var]])),
                     function(f){
   d <- valid_dat[get(facet_var)==f]
   #### Set min observation per group ####
   # d <- d[,x_var_n:=.N,by=x_var][x_var_n>1]
   #### Make celltypes an ordered factor ####
   if(is.numeric(d[[y_var]])){
     d[,median_y:=stats::median(get(y_var)),by=c(facet_var,x_var)]
     d[,(x_var):=factor(get(x_var),levels=unique(get(x_var)), ordered = TRUE)]
     data.table::setorderv(d, "median_y")
     method <- "ggbetweenstats"
   } else {
     d[,(y_var):=factor(get(y_var),levels=rev(unique(get(y_var))), ordered = TRUE)]
     counts <- d[,.N, by=c(x_var,y_var)][,pct:=.N/sum(N), by=c(x_var)][get(y_var)==levels(d[[y_var]])[1],]|>
       data.table::setorderv("pct",-1)
     d[,(x_var):=factor(get(x_var),levels=unique(counts[[x_var]]), ordered = TRUE)]
     method <- "ggbarstats"
   }

   tryCatch({
     if(method=="ggbetweenstats"){
       ggstatsplot::ggbetweenstats( data = d,
                                    x = !!ggplot2::sym(x_var),
                                    y = !!ggplot2::sym(y_var),
                                    title=f,
                                    type = type,
                                    centrality.plotting = centrality.plotting,
                                    centrality.type = centrality.type,
                                    pairwise.comparisons = pairwise.comparisons,
                                    xlab = x_var,
                                    ylab = y_var,
                                    ggtheme=
                                      ggstatsplot::theme_ggstatsplot() +
                                      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)),
                                    package = valid_pal$package[1],
                                    palette = valid_pal$palette[1],
                                    ...)
     } else {
       ggstatsplot::ggbarstats(data = d,
                               x = !!ggplot2::sym(y_var),
                               y = !!ggplot2::sym(x_var),
                               title=f,
                               type = type,
                               sample.size.label.args = list(size=3),
                               pairwise.comparisons = pairwise.comparisons,
                               xlab = y_var,
                               ylab = x_var,
                               ggtheme=
                                 ggstatsplot::theme_ggstatsplot() +
                                 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)),
                               ...)
     }
   }, error=function(e){message(e);NULL})
 })
 #### Filter out NULL ####
 plts_list <- plts_list[sapply(plts_list, function(x) !is.null(x))]
}

# #### Make celltypes an ordered factor ####
# plts <- ggstatsplot::grouped_ggbetweenstats( data = valid_dat,
#                                              x = !!ggplot2::sym(x_var),
#                                              y = !!ggplot2::sym(y_var),
#                                              grouping.var = !!ggplot2::sym(facet_var),
#                                              type = "nonparametric",
#                                              centrality.type = "nonparametric",
#                                              pairwise.comparisons = TRUE,
#                                              xlab = x_var,
#                                              ylab = y_var,
#                                              ggtheme=
#                                                ggstatsplot::theme_ggstatsplot() +
#                                                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)),
#                                              package = valid_pal$package[1],
#                                              palette = valid_pal$palette[1])

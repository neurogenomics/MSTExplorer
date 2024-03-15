#' Plot differential outcomes
#'
#' Plot differential outcomes (annotations) of different diseases as they occur
#' in different cell types via particular phenotypes.
#' @param remove_facets A character vector of facets to remove.
#' @param run_stats If \code{TRUE}, run statistical tests and display the
#' summary statistics directly on the plot.
#' @param max_facets Maximum number of facets to display.
#' @param q_threshold A list of thresholds for significance testing.
#'
#' @inherit plot_
#' @inheritParams plot_bar_dendro
#' @inheritParams KGExplorer::filter_dt
#' @inheritDotParams ggstatsplot::ggbetweenstats
#' @export
#' @examples
#' results <- add_symptom_results()
#' #### Multiple phenotypes per disease #####
#' results <- HPOExplorer::add_gpt_annotations(results)
#' p1 <- plot_differential_outcomes(results,
#'                            facet_var = "disease_name",
#'                            y_var = "severity_score_gpt")
plot_differential_outcomes <- function(results,
                                       filters = NULL,
                                       facet_var = "disease_name",
                                       x_var = "celltype_symptom",
                                       y_var = "severity_score_gpt",
                                       remove_facets = NULL,
                                       run_stats = FALSE,
                                       max_facets = NULL,
                                       q_threshold = list(
                                         summary = 0.05,
                                         pairwise = NULL
                                       ),
                                       ...
                                       ){
  requireNamespace("ggplot2")
  n_celltype <- n_ids <- q.value <- NULL;

  {
    results <- HPOExplorer::add_hpo_name(results)
    plot_dat <- results[!is.na(get(y_var)),]
    plot_dat[,n_celltype:=length(unique(get(x_var))),by=c(facet_var)]
    plot_dat[,n_ids:=.N, by=c(facet_var)]
    plot_dat <- plot_dat[n_celltype>1 & n_ids>1 & !is.na(get(x_var)),]
    plot_dat <- KGExplorer::filter_dt(dat = plot_dat,
                                      filters = filters)
  }

  if(isTRUE(run_stats)){
    #### Run pre-test to determine which models are sig ####
    # ggstatsplot::pairwise_comparisons()
    #### Run plot + stats ####
    plts <- plot_differential_outcomes_ggstatsplot(plot_dat=plot_dat,
                                                   remove_facets=remove_facets,
                                                   x_var=x_var,
                                                   y_var=y_var,
                                                   facet_var=facet_var,
                                                   max_facets=max_facets,
                                                   ...)
    stat_dat <- filter_ggstatsplot(plts=plts)
    stat_dat[["plot_data"]] <- plot_dat
    if(!is.null(q_threshold$summary)){
      n_sig <- nrow(stat_dat$summary_data[q.value<=q_threshold$summary])
      if(n_sig==0){
        warning("No significant results found. Returning unfiltered plot.")
        return(list(plot=plts,
                    data=stat_dat))
      }
    }
    #### Recreate the plot with only significant results ####
    if(!any(mapply(is.null,q_threshold))){
      plts2 <- plot_differential_outcomes_ggstatsplot(
        plot_dat=plot_dat[get(facet_var) %in% stat_dat$summary_data$facet,],
        remove_facets=remove_facets,
        x_var=x_var,
        y_var=y_var,
        facet_var=facet_var,
        max_facets=max_facets)
      stat_dat2 <- filter_ggstatsplot(plts=plts2,
                                      q_threshold=q_threshold)
      stat_dat2[["plot_data"]] <- plot_dat
      return(list(plot=plts2,
                  data=stat_dat2))
    } else {
      return(list(plot=plts,
                  data=stat_dat))
    }
  } else {
    requireNamespace("tidytext")

    if(!is.null(max_facets)){
      ids <- utils::head(unique(plot_dat[[facet_var]]),max_facets)
      plot_dat <- plot_dat[get(facet_var) %in% ids,]
    }
    #### Plot ####
    plt <- ggplot2::ggplot(plot_dat,
                           ggplot2::aes(x=tidytext::reorder_within(x=!!ggplot2::sym(x_var),
                                        by=!!ggplot2::sym(y_var),
                                        within=!!ggplot2::sym(facet_var)
           ),
           y=!!ggplot2::sym(y_var))) +
      ggplot2::geom_boxplot() +
      tidytext::scale_x_reordered() +
      ggplot2::labs(x=x_var) +
      ggplot2::facet_wrap(facets = facet_var,
                 nrow = 3,
                 scales = "free_x") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    return(list(plot=plt,
                data=plot_dat))
  }
}

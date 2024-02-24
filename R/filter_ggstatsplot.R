filter_ggstatsplot <- function(plts,
                               p_threshold = list(
                                 summary = NULL,
                                 pairwise = NULL
                               ),
                               q_threshold = list(
                                 summary = NULL,
                                 pairwise = NULL
                               )
                               ){

  stat_dat <- get_ggstatsplot_stats(plts = plts)
  summary_data <- stat_dat$summary_data
  pairwise_data <- stat_dat$pairwise_data
  p.value <- q.value <- NULL;
  #### Filter results ####
  if(!is.null(p_threshold$summary)){
    summary_data <- summary_data[p.value<=p_threshold$summary]
  }
  if(!is.null(q_threshold$summary)){
    summary_data <- summary_data[q.value<=q_threshold$summary]
  }
  if(!is.null(p_threshold$pairwise)){
    pairwise_data <- pairwise_data[p.value<=p_threshold$pairwise]
  }
  if(!is.null(q_threshold$pairwise)){
    pairwise_data <- pairwise_data[q.value<=q_threshold$pairwise]
  }
  return(list(summary_data=summary_data,
              pairwise_data=pairwise_data))
  # #### Remove specific plots from patchwork plts object ####
  # plots_to_remove <- stat_dat$summary_data$plot_id[
  #   !stat_dat$summary_data$plot_id %in% summary_data$plot_id
  # ]
  # keep <- summary_data$plot_id
  # ### Subset patchwork object
  # plts2 <- plts
  # for(i in plots_to_remove){
  #   plts2[i] = NULL
  # }
  # plts[[keep]]
  # ### remove specific subplots from patchwork object
  # patchwork::wrap_plots(plts)


}

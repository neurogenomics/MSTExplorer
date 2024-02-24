get_ggstatsplot_stats <- function(plts,
                                  method = "fdr"){
  plot_id <- NULL;
  #### patchwork object ####
  if(methods::is(plts,"patchwork") ||
     methods::is(plts,"list")){
    #### Get facet titles ####
    titles <- sapply(seq(length(plts)),
                     function(i){
                       plts[[i]]$labels$title
                     })
    #### Get model level-results ####
    summary_data <- lapply(seq(length(plts)),
                           function(i){
                             ggstatsplot::extract_stats(plts[[i]])$subtitle_data
                           }) |> data.table::rbindlist(idcol = "plot_id", fill=TRUE)
    summary_data[,facet:=titles[plot_id]]
    #### get pairwise comparisons results ####
    pairwise_data <- lapply(seq(length(plts)),
                            function(i){
                              ggstatsplot::extract_stats(plts[[i]])$pairwise_comparisons_data
                            }) |> data.table::rbindlist(idcol = "plot_id", fill=TRUE)
    pairwise_data[,facet:=titles[plot_id]]

    #### ggplot object ####
  } else if(methods::is(plts,"ggplot")){
    summary_data <- ggstatsplot::extract_stats(plts)$subtitle_data
    pairwise_data <- ggstatsplot::extract_stats(plts)$pairwise_comparisons_data
  } else {
    stopper("plts must be of class patchwork or ggplot")
  }
  #### Add multiple-testing correction ####
  summary_data$q.value <- stats::p.adjust(summary_data$p.value,
                                          method = method)
  pairwise_data$q.value <- stats::p.adjust(pairwise_data$p.value,
                                           method = method)
  #### Return ####
  return(list(summary_data = summary_data,
              pairwise_data = pairwise_data))
}

get_ggstatsplot_stats <- function(plts,
                                  method = "fdr"){
  plot_id <- facet <- NULL;
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
                              key_select <- c("pairwise_comparisons_data",
                                              "one_sample_data")
                              stats <- ggstatsplot::extract_stats(plts[[i]])
                              # Get the keys whose values in stats are not null and use the first one
                              {
                                key <- key_select[
                                  sapply(key_select, function(k){
                                    !is.null(stats[[k]])
                                  })][1]
                              }
                              stats[[key]]
                            }) |> data.table::rbindlist(idcol = "plot_id", fill=TRUE)
    pairwise_data[,facet:=titles[plot_id]]

    stats <- list(
      summary_data = summary_data,
      pairwise_data = pairwise_data
    )

    #### ggplot object ####
  } else if(methods::is(plts,"ggplot")){
    stats <- ggstatsplot::extract_stats(plts)
    summary_data <- data.table::data.table(stats$subtitle_data)
    # Get pairwise results
    if("pairwise_comparisons_data" %in% names(stats) &
       !is.null(stats$pairwise_comparisons_data) ){
      pairwise_data <- data.table::data.table(stats$pairwise_comparisons_data)
    } else if ("one_sample_data"  %in%  names(stats) &
               !is.null(stats$one_sample_data) ){
      pairwise_data <- data.table::data.table(stats$one_sample_data)
    } else{
      pairwise_data <- NULL
    }
  } else {
    stopper("plts must be of class patchwork or ggplot")
  }

  #### Add multiple-testing correction ####
  if(!is.null(summary_data)){
    summary_data$q.value <- stats::p.adjust(summary_data$p.value,
                                            method = method)
  }
  if(!is.null(pairwise_data)){
    pairwise_data$q.value <- stats::p.adjust(pairwise_data$p.value,
                                             method = method)
  }
  #### Return ####
  return(list(
    summary_data =  summary_data,
    pairwise_data =  pairwise_data
  ))
}

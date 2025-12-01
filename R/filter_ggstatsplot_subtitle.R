#' Filter ggstatsplot subtitle to show selected statistics only
#'
#' Filter ggstatsplot subtitle to show selected statistics only
#'
#' @param p ggstatsplot object
#' @param stats_idx Indices of statistics to keep in subtitle
#'
#' @return ggstatsplot object with filtered subtitle
#' @export
#' @examples
#'   p <- ggstatsplot::ggscatterstats(
#'     data = mtcars,
#'     x = disp,
#'     y = mpg
#'   )
#'   p_filtered <- filter_ggstatsplot_subtitle(p, stats_idx = c(1,3))
#'   print(p_filtered)
filter_ggstatsplot_subtitle <- function(p,
                                        stats_idx = c(1,3,4,7)){
    if (is.null(stats_idx)){
      return(p)
    }
    if (!1 %in% stats_idx){
        stop("stats_idx must include 1 to keep the list function.")
    }
    if (max(stats_idx) > length(p$labels$subtitle)){
        stop("stats_idx exceeds number of available statistics in subtitle.")
    }
    p$labels$subtitle <- p$labels$subtitle[stats_idx]
    return(p)
}

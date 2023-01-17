#' Count results
#'
#' Count the number of test that remain per celltype (or phenotype)
#' after filtering by one or more criterion.
#' @param group_var Variable to group results by.
#' @param sd_threshold The standard deviations from the mean
#' threshold to subset the results by.
#' @param p_threshold The p-value
#' threshold to subset the results by.
#' @inheritParams data.table::setorderv
#' @inheritParams ggnetwork_plot_full
#' @returns data.table with counts.
#'
#' @export
#' @importFrom data.table copy .N setorderv
#' @examples
#' results <- load_example_results()
#' celltype_counts <- count_results(results=results, group_var="CellType")
#' phenotype_counts <- count_results(results=results, group_var="Phenotype")
count_results <- function(results,
                          q_threshold = 0.0005,
                          fold_threshold = 1,
                          sd_threshold = NULL,
                          p_threshold = NULL,
                          group_var = "CellType",
                          order = -1){

  fold_change <- q <- p <- sd_from_mean <- NULL;

  dt <- data.table::copy(results)
  if(!is.null(q_threshold)){
    dt <- dt[q<=q_threshold,]
  }
  if(!is.null(fold_threshold)){
    dt <- dt[fold_change>=fold_threshold,]
  }
  if(!is.null(sd_threshold)){
    dt <- dt[sd_from_mean>=sd_threshold,]
  }
  if(!is.null(p_threshold)){
    dt <- dt[p>=p_threshold,]
  }
  counts <- dt[,.(count=.N), by=group_var]
  data.table::setorderv(counts,
                        cols = "count",
                        order = order)
  return(counts)
}

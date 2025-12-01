#' Calculate replication rate between two studies
#' @param dat A data.table containing p-values from two studies.
#' @param x The column name of the first study's p-values.
#' @param y The column name of the second study's p-values.
#' @param alpha Significance threshold for determining replication.
#' @return A list with replication rates from study A to B, B to A, and symmetric replication rate.
#' @examples
#' dat <- data.table::data.table(
#'   q_DescartesHuman = c(0.01, 0.2, 0.03, 0.5, 0.04),
#'   q_HumanCellLandscape = c(0.02, 0.03, 0.2, 0.6, 0.01)
#' )
#' replication_rates <- get_replication_rate(dat,
#' x="q_DescartesHuman", y="q_HumanCellLandscape",
#' )
#' print(replication_rates)
#' @export
get_replication_rate <- function(dat,
                                 x,
                                 y,
                                 alpha=0.05){
  sig_A <- sig_B <- NULL;

  dat[, sig_A := get(x) < alpha]
  dat[, sig_B := get(y) < alpha]

  RR_A_to_B <- dat[sig_A == TRUE, mean(sig_B)]
  RR_B_to_A <- dat[sig_B == TRUE, mean(sig_A)]
  RR_sym <- mean(c(RR_A_to_B, RR_B_to_A))

  return(
    list(
      RR_A_to_B = RR_A_to_B,
      RR_B_to_A = RR_B_to_A,
      RR_sym = RR_sym
    )
  )
}

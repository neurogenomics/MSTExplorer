#' Permutation test for correlation between two variables
#'
#' This function performs a permutation test to assess the significance of the correlation
#' between two variables in a data frame. It computes the observed correlation, generates
#' a null distribution of correlations by permuting one of the variables, and calculates
#' an empirical p-value.
#' @param df A data frame containing the variables of interest.
#' @param x_var A string specifying the name of the first variable (e.g., "info_content").
#' @param y_var A string specifying the name of the second variable (e.g., "n_cell_types_sig").
#' @param B An integer specifying the number of permutations to perform (default is 1000).
#' @param method A string specifying the correlation method to use ("pearson", "spearman", etc.; default is "pearson").
#' @param workers An integer specifying the number of parallel workers to use (default is 1).
#' @param seed An integer specifying the random seed for reproducibility (default is 1).
#' @param drop_duplicates A logical indicating whether to drop duplicate (x_var, y_var) pairs before analysis (default is TRUE).
#' @param save_path File path to cache results.
#' @param force_new Ignore cached results.
#' @inheritParams stats::cor
#' @return A data table containing the observed correlation, empirical p-value, and null distribution statistics.
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom stats cor quantile
#' @importFrom data.table data.table
#' @export
#' @examples
#' # Example usage:
#' df <- data.frame(
#'   info_content = rnorm(100),
#'   n_cell_types_sig = rnorm(100)
#' )
#' result <- run_permutation_test_cor(
#' df, x_var="info_content", y_var="n_cell_types_sig", B = 1000)
#' print(result)
run_permutation_test_cor <- function(df,
                                     x_var,      # e.g. "info_content"
                                     y_var,      # e.g. "n_cell_types_sig"
                                     B = 1000,
                                     method = "pearson",
                                     use = "complete.obs",
                                     workers = 1,
                                     seed = 1,
                                     drop_duplicates=FALSE,
                                     save_path=NULL,
                                     force_new=FALSE) {

  if (!is.null(save_path) && file.exists(save_path) && isFALSE(force_new)){
    messager("Loading cached results:",save_path)
    return(readRDS(save_path))
  }

  set.seed(seed)

  # Drop duplicates
  if (drop_duplicates){
    df <- df[!duplicated(df, by = c(x_var, y_var))]
  }

  # observed correlation (effect size)
  obs_res <- stats::cor.test(df[[x_var]], df[[y_var]],
                             method = method,
                             use = use)
  obs_cor <- obs_res$estimate
  obs_p <- obs_res$p.value

  # permutation distribution
  perm_cor <- BiocParallel::bplapply(
    seq_len(B),
    function(i) {
      perm_x <- sample(df[[x_var]])  # shuffle IC across phenotypes
      stats::cor(perm_x, df[[y_var]], method = method, use = "complete.obs")
    },
    BPPARAM = BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE)
  ) |> unlist()

  # two-sided empirical p-value
  p_emp <- mean(abs(perm_cor) >= abs(obs_cor))

  # optional: 95% null CI for the correlation
  null_ci <- stats::quantile(perm_cor, probs = c(0.025, 0.975), na.rm = TRUE)

  res <- data.table::data.table(
    x_var      = x_var,
    y_var      = y_var,
    method     = method,
    B          = B,
    obs_p      = obs_p,
    obs_cor    = obs_cor,
    p_emp      = p_emp,
    null_mean  = mean(perm_cor),
    null_ci_lo = null_ci[1],
    null_ci_hi = null_ci[2]
  )

  KGExplorer::cache_save(res, save_path = save_path)

  return(res)
}

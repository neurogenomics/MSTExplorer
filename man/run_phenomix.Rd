% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_phenomix.R
\name{run_phenomix}
\alias{run_phenomix}
\title{Run phenomix}
\usage{
run_phenomix(
  ctd_name,
  annotLevel = seq(length(ctd)),
  ymat,
  test_method = "glm",
  metric = "specificity",
  ctd = load_example_ctd(file = paste0("ctd_", ctd_name, ".rds")),
  save_path = file.path(tempfile(), "results", paste0("phenomix_", test_method, "_",
    metric), paste0("phenomix_", ctd_name, "_results.tsv.gz")),
  multivariate = FALSE,
  workers = NULL,
  force_new = FALSE,
  ...
)
}
\arguments{
\item{ctd_name}{Name of the CTD to load.}

\item{annotLevel}{An integer indicating which level of \code{sct_data} to
analyse (\emph{Default: 1}).}

\item{ymat}{gene x celltype matrix.}

\item{test_method}{Association testing method to use.}

\item{metric}{Which matrix within the CTD to use
(e.g. "mean_exp","specificity","specificity_quantiles").}

\item{ctd}{CellTypeDataset generated using
\link[EWCE]{generate_celltype_data}.}

\item{save_path}{Path to save the table of aggregated results to.}

\item{multivariate}{If \code{TRUE}, runs tests with each column in 
\code{xmat} as a multivariate predictor in a single model 
(one model per column in \code{yvar}). If\code{FALSE},
runs tests with each column in \code{xmat} as a univariate predictor in
separate models (one model per column in \code{yvar}).}

\item{workers}{Number (>1) or proportion (<1) of worker cores to use.}

\item{force_new}{Overwrite previous results
in the \code{save_dir_tmp}.}

\item{...}{Additional parameters passed to the statistical test function.}
}
\description{
Run many phenotype-cell type association tests in parallel using
\link[phenomix]{iterate_lm}.
}
\examples{
\dontrun{
  ymat <- HPOExplorer::hpo_to_matrix()
  ymat <- ymat[,1:10]
  lm_res <- run_phenomix(ctd_name = "HumanCellLandscape",
                         annotLevel = 3,
                         ymat = ymat,
                         save_path=tempfile())
}
}

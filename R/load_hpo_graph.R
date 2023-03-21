#' Load HPO graph
#'
#' Load the HPO as a precomputed graph.
#' @param file File to load.
#' @inheritParams load_example_ctd
#' @inheritParams piggyback::pb_download
#' @source
#' \code{
#' phenos <- HPOExplorer::make_phenos_dataframe(ancestor="All")
#' g <- HPOExplorer::make_igraph_object(phenos = phenos)
#' f <- file.path(tempdir(),"hpo_igraph.rds")
#' saveRDS(g,file = f)
#' piggyback::pb_upload(file = f,
#'                      tag = "v0.0.1", repo = "neurogenomics/MultiEWCE")
#' }
#' @returns graph object
#'
#' @keywords internal
#' @importFrom piggyback pb_download
#' @importFrom tools R_user_dir
#' @examples
#' g <- load_hpo_graph()
load_hpo_graph <- function(file="hpo_igraph.rds",
                           tag = "v0.0.1",
                           save_dir=tools::R_user_dir(package = "MultiEWCE")
                           ) {

  file <- file[[1]]
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  save_path <- file.path(save_dir,file)
  if (!file.exists(save_path)) {
    piggyback::pb_download(file = basename(file),
                           repo = "neurogenomics/MultiEWCE",
                           tag = tag,
                           dest = save_dir,
                           overwrite = TRUE)
  }
  g <- readRDS(save_path)
  return(g)
}

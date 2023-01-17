#' Load example CTD (Descartes)
#'
#' This loads a example CTD for vignettes and testing purposes.
#' It is \href{https://descartes.brotmanbaty.org/}{Descartes} Gene Expression
#' During Development single-cell RNA-seq atlas at annotation level 1.
#' @param file Name of a remotely stored file.
#' @param save_dir Where to store the file locally.
#' @inheritParams piggyback::pb_download
#'
#' @export
#' @importFrom piggyback pb_download
#' @importFrom tools R_user_dir
#' @examples
#' CTD <- load_example_ctd()
load_example_ctd <- function(file="CTD_Descartes_example.rds",
                             tag = "v0.0.1",
                             save_dir=tools::R_user_dir(package = "MultiEWCE")
                             ) {

  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  save_path <- file.path(save_dir,file)
  if (!file.exists(save_path)) {
    piggyback::pb_download(file = basename(file),
                           repo = "neurogenomics/MultiEWCE",
                           tag = tag,
                           dest = save_dir,
                           overwrite = TRUE)
  }
  return(readRDS(save_path))
}

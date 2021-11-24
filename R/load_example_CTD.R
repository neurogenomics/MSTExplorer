#' Load example CTD (Descartes)
#'
#' This loads a example CTD for vignettes and testing purposes. It is Descartes
#' human cell atlas data at annotation level 1.
#'
#' @examples
#' CTD <- load_example_CTD()
#'
#' @export
load_example_CTD <- function() {
  if (! file.exists("CTD_Descartes_example.rds")) {
    piggyback::pb_download(repo="neurogenomics/MultiEWCE",
                           tag = "v0.0.1")

    return(readRDS("CTD_Descartes_example.rds"))
  } else {
    return(readRDS("CTD_Descartes_example.rds"))
  }
}

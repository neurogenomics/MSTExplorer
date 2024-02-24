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
#' @examples
#' CTD <- load_example_ctd()
load_example_ctd <- function(file= c("ctd_DescartesHuman_example.rds",
                                     "ctd_DescartesHuman.rds",
                                     "ctd_HumanCellLandscape.rds"
                                     ),
                             multi_dataset = FALSE,
                             tag = "latest",
                             save_dir=KGExplorer::cache_dir(package="MSTExplorer")
                             ) {

  if(isFALSE(multi_dataset)){
    file <- file[1]
  }
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  ctd_list <- lapply(stats::setNames(file,
                                     gsub("\\.rds$|ctd_","",file)),
                     function(f){
    messager("Loading",f)
    save_path <- file.path(save_dir,f)
    if (!file.exists(save_path)) {
      piggyback::pb_download(file = basename(f),
                             repo = "neurogenomics/MSTExplorer",
                             tag = tag,
                             dest = save_dir,
                             overwrite = TRUE)
    }
    return(readRDS(save_path))
  })
  if(isFALSE(multi_dataset)){
    ctd_list <- ctd_list[[1]]
  }
  return(ctd_list)
}

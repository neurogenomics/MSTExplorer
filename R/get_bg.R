#' Get background genes
#'
#' Generate background genes given one or more species.
#' Caches the list to avoid excessive API calls to
#' \href{https://biit.cs.ut.ee/gprofiler}{g:Profiler}.
#' @inheritParams get_data
#' @inheritParams orthogene::create_background
#' @inheritDotParams orthogene::create_background
#' @returns A vector of background genes.
#'
#' @export
#' @examples
#' bg <- get_bg()
get_bg <- function(species1 = "human",
                   species2 = "human",
                   method = "gprofiler",
                   save_dir = KGExplorer::cache_dir(package="MultiEWCE"),
                   overwrite = FALSE,
                   verbose = TRUE,
                   ...){
  requireNamespace("orthogene")
  #### Create save path ####
  save_path <- file.path(
    save_dir,
    paste0(paste("bg",species1,species2,method,sep="-"),".rds")
    )
  if(file.exists(save_path) &&
     isFALSE(overwrite)){
    messager("Useing cached bg.",v=verbose)
    bg <- readRDS(save_path)
  } else {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    bg <- orthogene::create_background(species1 = species1,
                                       species2 = species2,
                                       method = method,
                                       ...)
    ## Add version attribute using date
    attr(bg,"version") <- as.character(Sys.Date())
    saveRDS(bg,save_path)
  }
  get_version(obj = bg,
              verbose = verbose)
  return(bg)
}

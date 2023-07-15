#' Get remote data
#'
#' Download remotely stored data via \link[piggyback]{pb_download}.
#' @keywords internal
#' @inheritParams piggyback::pb_download
#' @importFrom tools R_user_dir
get_data <- function(fname,
                     repo = "neurogenomics/MultiEWCE",
                     storage_dir = tools::R_user_dir(
                       package = "MultiEWCE",
                       which = "cache"
                     ),
                     overwrite = TRUE,
                     tag = "latest",
                     check = FALSE
                     ){
  requireNamespace("piggyback")
  Sys.setenv("piggyback_cache_duration" = 10)

  tmp <- file.path(storage_dir, fname)
  dir.create(storage_dir, showWarnings = FALSE, recursive = TRUE)
  piggyback::pb_download(
    file = fname,
    dest = storage_dir,
    repo = repo,
    overwrite = overwrite,
    tag = tag
  )
  #### Read/return ####
  if(grepl("\\.rds$",tmp)){
    return(readRDS(tmp))
  } else if(grepl("\\.csv$|\\.tsv$",tmp)){
    return(data.table::fread(tmp))
  } else {
    return(tmp)
  }
}

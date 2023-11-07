#' Get remote data
#'
#' Download remotely stored data via \link[piggyback]{pb_download}.
#' @param save_dir Directory to save data to.
#' @keywords internal
#' @inheritParams piggyback::pb_download
#' @importFrom tools R_user_dir
get_data <- function(fname,
                     repo = "neurogenomics/MultiEWCE",
                     save_dir = tools::R_user_dir(
                       package = "MultiEWCE",
                       which = "cache"
                     ),
                     overwrite = TRUE,
                     tag = "latest",
                     check = FALSE
                     ){
  requireNamespace("piggyback")
  Sys.setenv("piggyback_cache_duration" = 10)

  tmp <- file.path(save_dir, fname)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  piggyback::pb_download(
    file = fname,
    dest = save_dir,
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

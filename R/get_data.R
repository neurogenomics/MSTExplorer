#' Get remote data
#'
#' Download remotely stored data via \link[piggyback]{pb_download}.
#' @param save_dir Directory to save data to.
#' @keywords internal
#' @inheritParams piggyback::pb_download
#' @inheritDotParams piggyback::pb_download
get_data <- function(fname,
                     repo = "neurogenomics/MSTExplorer",
                     save_dir = KGExplorer::cache_dir(package="MSTExplorer"),
                     overwrite = TRUE,
                     tag = "latest",
                     check = FALSE,
                     ...
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
    tag = tag,
    ...
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

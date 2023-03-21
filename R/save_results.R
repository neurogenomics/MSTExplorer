save_results <- function(results,
                         save_dir,
                         prefix="results_",
                         verbose=TRUE){
  if (!is.null(save_dir) &&
      nrow(results)>0) {
    save_path <- file.path(
      save_dir,
      gsub(" ","_",
           paste0(prefix,stringr::str_replace_all(Sys.time(),":","-"),
                  ".rds")
      )
    )
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    messager("\nSaving results ==>",save_path,v=verbose)
    saveRDS(results,save_path)
    return(save_path)
  }
}

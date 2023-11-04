save_results <- function(results,
                         save_path,
                         verbose=TRUE){

  if (!is.null(save_path) &&
      nrow(results)>0) {
    messager("\nSaving results ==>",save_path,v=verbose)
    saveRDS(results,save_path)
    return(save_path)
  }
}

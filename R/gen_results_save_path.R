gen_results_save_path <- function(save_dir,
                                  prefix="gen_results",
                                  suffix=NULL){
  save_path <- file.path(
    save_dir,
    gsub(" ","_",paste0(prefix,suffix,".rds"))
  )
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  return(save_path)
}

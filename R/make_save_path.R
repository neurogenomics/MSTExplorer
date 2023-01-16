make_save_path <- function(save_dir,
                           list_name){
  file.path(save_dir,
            paste0(gsub(" +","_",tolower(list_name)),".rds"))
}

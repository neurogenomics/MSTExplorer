#' Merge results
#'
#' Merges individual .rds files from the results output directory into one.
#' dataframe. Can also handle results lists directly.
#' Alternatively, can simply provide \code{save_dir} and \emph{.rds}
#' files will be searched for and imported.
#' @param res_files The output of \link[MultiEWCE]{ewce_para}
#' (can be either a list of file names, or a nested list of results).
#' @inheritParams ewce_para
#' @inheritParams gen_results
#' @returns dataframe and merged results.
#'
#' @export
#' @importFrom data.table rbindlist
#' @examples
#' gene_data <- HPOExplorer::load_phenotype_to_genes()
#' ctd <- MultiEWCE::load_example_ctd()
#' list_names <- unique(gene_data$hpo_name)[seq(3)]
#' res_files <- ewce_para(ctd = ctd,
#'                        gene_data = gene_data,
#'                        list_names = list_names,
#'                        reps = 10)
#' all_results <- merge_results(res_files=res_files)
merge_results <- function(save_dir=NULL,
                          res_files=NULL,
                          list_name_column = "hpo_name") {

  if(is.null(res_files)){
    if(is.null(save_dir)) stop("Must provided save_dir when res_files=NULL.")
    res_files <- list.files(path = save_dir,
                            pattern = ".rds$",
                            ignore.case = TRUE,
                            full.names = TRUE)
    names(res_files) <- gsub("_"," ",tolower(gsub(".rds$","",res_files)))
    messager(formatC(length(res_files),big.mark = ","),"results files found.")
  }
  lapply(seq(length(res_files)),
         function(i){
   if(is.null(res_files[[i]])){
     return(NULL)
   } else if(is.list(res_files[[i]])){
     cur <- res_files[[i]]$results
   } else{
     cur <- readRDS(res_files[[i]])$results
   }
    cur[[list_name_column]] <- names(res_files)[[i]]
    return(cur)
  }) |> data.table::rbindlist(fill=TRUE)
}

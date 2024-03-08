#' Run phenomix
#'
#'
#' @export
#' @examples
#' ymat <- HPOExplorer::hpo_to_matrix()
#' ymat <- ymat[,1:10]
#' lm_res <- run_phenomix(ctd_name = "HumanCellLandscape",
#'                        annotLevel = 3,
#'                        ymat = ymat,
#'                        save_path=tempfile())
run_phenomix <- function(ctd_name,
                         annotLevel,
                         ymat,
                         test_method="glm",
                         metric="specificity",
                         ctd = MSTExplorer::load_example_ctd(
                           file = paste0("ctd_",ctd_name,".rds")
                           ),
                         xmat = ctd[[annotLevel]][[metric]],
                         save_path = here::here(
                           "results",paste0("phenomix_",test_method,"_",metric),
                           paste0("phenomix_",ctd_name,"_results.tsv.gz")
                         ),
                         multivariate=FALSE,
                         workers = NULL,
                         force_new = FALSE,
                         ...
){
  if(file.exists(save_path) && isFALSE(force_new)){
    message("Loading existing results from ",save_path)
    return(data.table::fread(save_path))
  }
  if(test_method=="glm_univariate"){
    test_method <- "glm"
    multivariate <- FALSE
  }
  lm_res <- phenomix::iterate_lm(xmat = xmat,
                                 ymat = ymat,
                                 test_method = test_method,
                                 multivariate = multivariate,
                                 workers = workers,
                                 ...)
  run_phenomix_postprocess <- function(lm_res,
                                       annotLevel,
                                       effect_var = c("estimate","statistic",
                                                      "ges","F")){
    effect_var <- intersect(effect_var,names(lm_res))[1]
    lm_res[,annotLevel:=annotLevel]
    data.table::setnames(lm_res,
                         c("xvar","yvar"),
                         c("CellType","hpo_id"))
    # x-u/sd
    # lm_res[,tmp:=log(scales::rescale(estimate,c(1,2)))][,fold_change:=(tmp-mean(tmp))/sd(tmp)]
    # lm_res[,fold_change:=log(abs(estimate)/mean(abs(estimate)))][,fold_change:=fold_change+abs(min(fold_change))]
    # lm_res[,fold_change:=scales::rescale(log(abs(estimate)),c(0,5))]
    lm_res[,fold_change:=log(scales::rescale(get(effect_var),c(1,5)))]
    # hist(lm_res$fold_change)
  }
  run_phenomix_postprocess(lm_res,
                           annotLevel=annotLevel)
  #### Save ####
  if(!is.null(save_path)){
    messager("Saving results -->",save_path)
    dir.create(dirname(save_path), recursive = TRUE,showWarnings = FALSE)
    data.table::fwrite(lm_res,save_path)
  }
  return(lm_res)
}




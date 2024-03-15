#' Run phenomix
#'
#' Run many phenotype-cell type association tests in parallel using
#' \link[phenomix]{iterate_lm}.
#' @param ctd_name Name of the CTD to load.
#' @param metric Which matrix within the CTD to use
#' (e.g. "mean_exp","specificity","specificity_quantiles").
#' @param save_path Path to save the table of aggregated results to.
#' @inheritParams ewce_para
#' @inheritParams phenomix::iterate_lm
#' @export
#' @examples
#' \dontrun{
#'   ymat <- HPOExplorer::hpo_to_matrix()
#'   ymat <- ymat[,1:10]
#'   lm_res <- run_phenomix(ctd_name = "HumanCellLandscape",
#'                          annotLevel = 3,
#'                          ymat = ymat,
#'                          save_path=tempfile())
#' }
run_phenomix <- function(ctd_name,
                         annotLevel,
                         ymat,
                         test_method="glm",
                         metric="specificity",
                         ctd = load_example_ctd(
                           file = paste0("ctd_",ctd_name,".rds")
                           ),
                         xmat = ctd[[annotLevel]][[metric]],
                         save_path = file.path(
                           tempfile(),
                           # here::here(),
                           "results",paste0("phenomix_",test_method,"_",metric),
                           paste0("phenomix_",ctd_name,"_results.tsv.gz")
                         ),
                         multivariate=FALSE,
                         workers = NULL,
                         force_new = FALSE,
                         ...
){
  requireNamespace("phenomix")
  effect <- NULL;

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
    # lm_res[,tmp:=log(scales::rescale(estimate,c(1,2)))][,effect:=(tmp-mean(tmp))/stats::sd(tmp)]
    # lm_res[,effect:=log(abs(estimate)/mean(abs(estimate)))][,effect:=effect+abs(min(effect))]
    # lm_res[,effect:=scales::rescale(log(abs(estimate)),c(0,5))]
    lm_res[,effect:=log(scales::rescale(get(effect_var),c(1,5)))]
    # hist(lm_res$effect)
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




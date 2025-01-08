#' Test target cell types
#'
#' Test the following hypotheses across \link{gen_results} enrichment results.
#' @param tests The types of tests to conduct.
#' \itemize{
#' \item{"within_branches"}{
#' Are on-target cell types enriched for significant results?}
#' \item{"across_branches"}{
#' Are on-target cell types enriched for significant results across branches?}
#' \item{"across_branches_per_celltype"}{
#' Are on-target cell types enriched for significant results across
#' branches per cell type?}
#' }
#' @param ancestor_var The name of the \code{results} column
#'  containing the ancestor name.
#' @param within_var Within-subject variable.
#' @param method Methods to run stats with.
#' @inheritParams ggnetwork_plot_full
#' @inheritParams plot_bar_dendro
#' @inheritParams plot_bar_dendro_facets
#' @inheritParams ewce_para
#' @inheritParams rstatix::anova_test
#' @export
#' @examples
#' res <- test_target_celltypes(tests="within_branches")
test_target_celltypes <- function(results=load_example_results(),
                                  target_celltypes = get_target_celltypes(),
                                  celltype_col="cl_id",
                                  tests=c("within_branches",
                                          "across_branches",
                                          "across_branches_per_celltype"),
                                  method=c("glm","anova"),
                                  within_var="hpo_id",
                                  ancestor_var="ancestor_name",
                                  q_threshold=0.05,
                                  cores=1
                                  ){
  requireNamespace("dplyr")
  requireNamespace("rstatix")
  ancestor_name <- is_sig <- is_target <- valid <- cl_id <- NULL;

  method <- match.arg(method)
  results <- HPOExplorer::add_ancestor(results)
  results <- map_celltype(results)
  results[,is_sig:=q<q_threshold][,is_target:=get(celltype_col) %in%
                                  target_celltypes[[get(ancestor_var)]],
                                  by=ancestor_var]

  res <- list()
  #### Are on-target cell types enriched for significant results?
  if("within_branches" %in% tests){
    messager("Running tests: within_branches")
    res[["within_branches"]] <- results|>
      dplyr::filter(!!dplyr::sym(ancestor_var) %in% names(target_celltypes))|>
      dplyr::group_by(!!dplyr::sym(ancestor_var))|>
      rstatix::anova_test(formula = is_sig ~ is_target,
                          # between = "is_target",
                          within=dplyr::all_of(within_var)
                          ) |>
      dplyr::ungroup() |>
      rstatix::adjust_pvalue(method = "bonf") |>
      rstatix::add_significance()
  }


  #### Are on-target cell types enriched for significant results across branches?
  if("across_branches" %in% tests){
    messager("Running tests: across_branches")
    res[["across_branches"]] <- parallel::mclapply(
      stats::setNames(names(target_celltypes),
                      names(target_celltypes)),
      function(b){
        messager("Running:",b, parallel = TRUE)
        d <- data.table::copy(results)
        d[,is_sig:=q<0.05][,is_target:=get(celltype_col) %in%
                             target_celltypes[b][[ancestor_name]], by=.I]
        d[get(celltype_col) %in% unlist(target_celltypes),]|>
          # dplyr::group_by(cl_name)|>
          rstatix::anova_test(formula = is_sig ~ is_target,
                              within=dplyr::all_of(within_var)
                              )
      }, mc.cores = cores) |> data.table::rbindlist(idcol = "ancestor_name")|>
      rstatix::adjust_pvalue(method = "bonf") |>
      rstatix::add_significance()
  }


  #### For branch per celltype, are the on-target celltypes more often enriched
  ## in the respective on-target branch than any other branch?
  if("across_branches_per_celltype" %in% tests){
    messager("Running tests: across_branches_per_celltype")
    res[["across_branches_per_celltype"]] <- parallel::mclapply(
      stats::setNames(names(target_celltypes),
                      names(target_celltypes)),
      function(b){
        messager("Running tests: ",b, parallel = TRUE)
        d <- data.table::copy(results)
        d[,is_sig:=q<0.05]
        ## Define on-target cell types for branch
        d[,is_target:=(get(celltype_col) %in% target_celltypes[[b]]) &
            (ancestor_name==b)]
        ### Ensure there's enough variation to run tests
        d[,valid:=(length(unique(is_target))>1) & (length(unique(is_sig))>1),
          by=c(celltype_col)]
        d <- d[valid==TRUE]
        if(nrow(d)==0){
          messager("Skipping tests.")
          return(NULL)
        }
        d[,is_sig:=q<0.05][,is_target:=get(celltype_col) %in%
                             target_celltypes[b][[ancestor_name]], by=.I]
        if(method=="anova"){
          d2 <- d|>
            dplyr::group_by(cl_id)|>
            rstatix::anova_test(formula = is_sig ~ is_target,
                                within=dplyr::all_of(within_var)
            )
        } else{
          ## This method becomes horribly slow (never finishes)
          ## if you try to include hpo_id as a within-subject variable.
          lapply(unique(d[[celltype_col]]), function(x){
            d2 <- d[get(celltype_col)==x]|>
              stats::glm(formula= is_sig ~ is_target,
                         family="binomial") |>
              rstatix::tidy()
            d2[[celltype_col]] <- x
            d2
          })|>data.table::rbindlist(fill=TRUE)
        }
        # lme4::lmer(data=d[get(celltype_col)==x],
        #            formula= is_sig ~ is_target + (1|hpo_id),
        #            control=lme4::lmerControl(check.nobs.vs.nlev = "ignore",
        #                                      check.nobs.vs.rankZ = "ignore",
        #                                      check.nobs.vs.nRE="ignore")) |>
        #   summary() |>as.data.frame()
      }, mc.cores = cores) |>
      data.table::rbindlist(idcol = "ancestor_name")|>
      rstatix::adjust_pvalue(method = "bonf") |>
      rstatix::add_significance()
  }
  #### Return ####
  return(res)
}

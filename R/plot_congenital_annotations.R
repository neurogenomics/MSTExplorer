#' Plot congenital annotations
#'
#' Test whether there is a difference in proportion if significantly associated
#'  fetal and non-fetal cell types across phenotypes with congenital onset
#'  vs. those without.
#' @inheritParams prioritise_targets
#' @param gpt_annot A data.table of GPT annotations.
#' @param fetal_keywords A character vector of keywords to identify fetal cell types.
#' @param celltype_col The column name of the cell type.
#' @param remove_annotations A character vector of annotations to remove.
#' @param save_path The path to save the plot.
#' @param by_branch Use HPO ancestors as the x-axis instead of the
#'  frequency of congenital onset.
#' @param x_var X-axis variable to plot.
#' @param add_baseline Add a horizontal line showing the proportions
#'  expected by random.
#' @inheritParams ggstatsplot::ggbarstats
#' @inheritDotParams ggstatsplot::ggbarstats
#' @export
#' @examples
#' results <- load_example_results()
#' results2 <- plot_congenital_annotations(results=results)
plot_congenital_annotations <- function(results,
                                        gpt_annot = HPOExplorer::gpt_annot_codify(),
                                        hpo=HPOExplorer::get_hpo(),
                                        # fetal_keywords=NULL,
                                        fetal_keywords=c("fetal",
                                                         "fetus",
                                                         "primordial",
                                                         "hESC",
                                                         "embryonic"),
                                        celltype_col="author_celltype",
                                        x_var=c("fetal_celltype",
                                                "fetal_only"),
                                        remove_annotations=c("varies"),
                                        keep_descendants=NULL,
                                        by_branch=FALSE,
                                        q_threshold=0.05,
                                        package="palettetown",
                                        palette="mewtwo",
                                        proportion.test	= TRUE,
                                        add_baseline=TRUE,
                                        save_path=NULL,
                                        ...){
  requireNamespace("ggstatsplot")
  fetal_celltype <- ancestor_name <- congenital_onset <- has_adult_and_fetal <-
    NULL;
  x_var <- match.arg(x_var)
  {
    results <- prepare_congenital_annotations(results = results,
                                              fetal_keywords = fetal_keywords,
                                              celltype_col = celltype_col,
                                              gpt_annot = gpt_annot)
    #### Filter ####
    dat <- results[q<q_threshold & (!congenital_onset %in% remove_annotations),]
    #### Add fetal_only col ####
    dat[,fetal_only:=(all(stage!="Adult") & has_adult_and_fetal==TRUE),
        by=c("hpo_id","cl_name")][has_adult_and_fetal==FALSE,fetal_only:=NA]
    #### Add ancestors ####
    dat <- HPOExplorer::add_ancestor(dat,
                                     keep_descendants = keep_descendants,
                                     hpo = hpo)
  }
  #### Bar plot ####
  if(by_branch){
    #### Branches plot ####
    dat[,c("fetal_celltype_prop","fetal_only_prop"):=list(
      mean(fetal_celltype),
      mean(fetal_only)
    ),by=c("ancestor_name")]
    data.table::setorderv(dat,"fetal_celltype_prop",-1)
    dat[,ancestor_name:=factor(ancestor_name,
                               levels=unique(ancestor_name),
                               ordered = TRUE)]
    if(x_var=="fetal_only"){
      ggbar <- ggstatsplot::ggbarstats(
        dat,
        x="fetal_only",
        y="ancestor_name",
        package=package,
        palette=palette,
        proportion.test=proportion.test,
        ...)
    } else {
      ggbar <- ggstatsplot::ggbarstats(
        dat,
        x="fetal_celltype",
        y="ancestor_name",
        package=package,
        palette=palette,
        proportion.test=proportion.test,
        ...)
    }
    ggbar <- ggbar +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 1, vjust=.5))
  } else{
    if(x_var=="fetal_only"){
      ggbar <- ggstatsplot::ggbarstats(
        dat,
        x="fetal_only",
        y="congenital_onset",
        package=package,
        palette=palette,
        proportion.test=proportion.test,
        ...)
    } else {
      ggbar <- ggstatsplot::ggbarstats(
        dat,
        x="fetal_celltype",
        y="congenital_onset",
        package=package,
        palette=palette,
        proportion.test=proportion.test,
        ...)
    }
  }
  #### Add baseline ####
  if(isTRUE(add_baseline)){
    if(x_var %in% names(results)){
      baseline <- (table(results[[x_var]])/nrow(results))[1]
    } else {
      baseline <- (table(dat[[x_var]])/nrow(dat))[1]
    }
     ggbar <- ggbar +
       ggplot2::geom_hline(yintercept=baseline,
                           linetype="dashed",
                           color="black", alpha=.8) +
       ggplot2::annotate("text",
                         size=3,
                         x=1,
                         y=baseline,
                         label=paste0("Baseline: ",round(baseline*100),"%"),
                         hjust=0.8,
                         vjust=-1)
   }
  ## Extract test results
  data_stats <- get_ggstatsplot_stats(ggbar)
  #### Save ####
  KGExplorer::plot_save(ggbar, save_path)
  #### Return ####
  return(
    list(
      raw_data=results,
      data=ggbar$data,
      data_stats=data_stats,
      plot=ggbar
    )
  )
}




# Attempt to make pvalues not exactly 0
# format_value_og <- insight::format_value
# format_value <- function(value,digits,...){
#   sapply(value,function(x){
#     if(x==0){
#       x <- .Machine$double.xmin
#     }
#     format_value_og(x, digits=digits,...)
#   })
# }
# assignInNamespace("format_value",format_value, "insight")
# #### Replace exact 0s ####
# ggb <- ggplot2::ggplot_build(ggbar)
# terms(eval(rlang::parse_exprs(ggb$data[[3]]$label[1])[[1]])[[1]])
#
# ggb$plot$layers[[4]]$data <-
#   ggb$plot$layers[[4]]$data|>
#   dplyr::filter(p.value==0)|>
#   dplyr::mutate(p.label=glue::glue("list(~italic(p)=='{insight::format_value(.Machine$double.xmin, ggbar$plot_env$digits)}')"))
# ggplot2::ggplot_gtable(ggb)


## Test relationship with severity_class
# gpt_class <- HPOExplorer::gpt_annot_class()
# dat2 <- data.table::merge.data.table(
#   dat,
#   gpt_class[,c("hpo_id","severity_score_gpt","severity_class")],
#   allow.cartesian = TRUE)
# ggbar <- ggstatsplot::ggbarstats(
#   dat2,
#   x="fetal_only",
#   y="severity_class",
#   package=package,
#   palette=palette,
#   proportion.test=proportion.test,
#   ...)

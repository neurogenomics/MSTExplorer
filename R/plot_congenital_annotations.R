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
#' @export
#' @examples
#' results <- load_example_results()
#' results2 <- add_ctd(results=results)
plot_congenital_annotations <- function(results,
                                        gpt_annot = HPOExplorer::gpt_annot_codify(),
                                        fetal_keywords=c("fetal",
                                                         "fetus",
                                                         "primordial",
                                                         "hESC",
                                                         "embryonic"),
                                        celltype_col="author_celltype",
                                        remove_annotations=c("varies"),
                                        save_path=NULL){
  requireNamespace("ggstatsplot")
  fetal_celltype <- congenital_onset <- NULL;

  results <- HPOExplorer::add_gpt_annotations(results,
                                              annot = gpt_annot$annot)
  results <- map_celltype(results)
  results[,fetal_celltype:=grepl(paste(fetal_keywords,collapse="|"),
                                 get(celltype_col),
                                 ignore.case = TRUE)]
  ggbar <- ggstatsplot::ggbarstats(
    results[q<0.05 & (!congenital_onset %in% remove_annotations),],
    x="fetal_celltype",
    y="congenital_onset",
    package="palettetown",
    palette="mewtwo",
    proportion.test	= TRUE)
  #### Extract test results ####
  data_stats <- get_ggstatsplot_stats(ggbar)
  #### Save ####
  KGExplorer::plot_save(ggbar, save_path)
  #### Return ####
  return(
    list(
      data=ggbar$data,
      data_stats=data_stats,
      plot=ggbar
    )
  )
}

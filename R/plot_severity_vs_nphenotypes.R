#' Plot severity vs. number of phenotypes
#'
#' Plot the mean composite GPT severity score of all all the phenotypes each
#' cell type is significantly associated with vs. the number of phenotypes
#' each cell types is significantly associated with.
#'
#' @inheritParams prioritise_targets
#' @inheritParams ggrepel::geom_label_repel
#' @inheritDotParams ggstatsplot::ggscatterstats
#' @param gpt_annot A data.table of GPT annotations.
#' @param size Label text size.
#' @param n_label The number of top and bottom cell types to label in the plot.
#' Top/bottom terms are determined by sorting on both the x and y axes.
#' @param remove_subtitle Remove the formula in the subtitle.
#' See \link[ggstatsplot]{ggscatterstats} for details.
#' @param remove_caption Remove the formula in the caption.
#' See \link[ggstatsplot]{ggscatterstats} for details.
#' @inheritParams ggplot2::theme_bw
#' @export
#' @examples
#' results <- load_example_results()
#' out <- plot_severity_vs_nphenotypes(results=results)
plot_severity_vs_nphenotypes <- function(results,
                                         gpt_annot=HPOExplorer::gpt_annot_codify(),
                                         cl = KGExplorer::get_ontology("cl", remove_rings=TRUE),
                                         q_threshold=0.5,
                                         n_label=3,
                                         size=3,
                                         min.segment.length = 0,
                                         remove_subtitle=FALSE,
                                         remove_caption=TRUE,
                                         point_fill=ggplot2::alpha("white",.75),
                                         base_size=8,
                                         run_prune_ancestors=FALSE,
                                         ...){
  requireNamespace("ggstatsplot")
  severity_score_gpt <- hpo_id <- phenotypes_per_celltype <- .I <- NULL;

  ## Merge and annotate results
  severity_vs_npheno <- merge(gpt_annot$annot_weighted,
                              map_celltype(results[q<q_threshold])
                              )[,list(
        severity_score_gpt=mean(severity_score_gpt, na.rm=TRUE),
        phenotypes_per_celltype=data.table::uniqueN(hpo_id),
        "log10(severity_score_gpt)"=log10(mean(severity_score_gpt, na.rm=TRUE)),
        "log10(phenotypes_per_celltype)"=log10(data.table::uniqueN(hpo_id))
                              ),
                              by=c("cl_name","cl_id")]
  # Prune redundant cell types
  if(run_prune_ancestors){
    severity_vs_npheno <- KGExplorer::prune_ancestors(dat = severity_vs_npheno,
                                                      id_col = "cl_id",
                                                      ont = cl)
  }
  # Plot
  severity_vs_npheno_plot <- ggstatsplot::ggscatterstats(
    severity_vs_npheno,
     x="log10(severity_score_gpt)",
     y="log10(phenotypes_per_celltype)",
     xsidehistogram.args = list(color=ggplot2::alpha("black",.5),
                                fill=ggplot2::alpha("#8f48e7",.5)),
     ysidehistogram.args = list(color=ggplot2::alpha("black",.5),
                                fill=ggplot2::alpha("#8f48e7",.5)),
    ...
  ) +
    ggplot2::labs(x=expression(log[10]~"(GPT severity score)"),
                  y=expression(log[10]~"(Phenotypes per celltype)")) +
    ggplot2::theme_bw(base_size = base_size)
  if(remove_subtitle){
    severity_vs_npheno_plot$labels$subtitle <- NULL
  }
  if(remove_caption){
    severity_vs_npheno_plot$labels$caption <- NULL
  }
  # severity_vs_npheno_plot$labels$group <- NULL
  ## Top/bottom celltype labels
  severity_vs_npheno_plot_data <- data.table::copy(severity_vs_npheno_plot$data)
  severity_vs_npheno_plot_data[,xy_mean:=mean(
    c(log10(severity_score_gpt),
      log10(phenotypes_per_celltype))),
                               by=.I] |>
    data.table::setorderv("xy_mean",-1)
  top_bottom_dat <- rbind(
    severity_vs_npheno_plot_data |>
      dplyr::slice_max(dplyr::tibble(severity_score_gpt,
                                     phenotypes_per_celltype), n = n_label),
    severity_vs_npheno_plot_data |>
      dplyr::slice_min(dplyr::tibble(severity_score_gpt,
                                     phenotypes_per_celltype), n = n_label)
  )
  gg <- severity_vs_npheno_plot +
    ggrepel::geom_label_repel(data=top_bottom_dat,
                              min.segment.length = min.segment.length,
                              size=size,
                              fill=point_fill,
                              ggplot2::aes(label=cl_name))
  ## Return
  return(
    list(data=severity_vs_npheno,
         plot=gg)
  )
}

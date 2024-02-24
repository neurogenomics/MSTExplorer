#' Plot bar dendrogram
#'
#' Create a plot summarising MSTExplorer results as a bar chart with multiple
#'  facets and a cell ontology-based dendrogram.
#' @param celltype_col Name of the cell type column in the \code{results}.
#' @inheritParams ggnetwork_plot_full
#' @inheritParams HPOExplorer::add_ancestor
#' @inheritParams ggplot2::facet_wrap
#' @inheritParams ggplot2::theme
#' @inheritDotParams EWCE::ewce_plot
#' @returns A bar chart with dendrogram of EWCE results in each cell type.
#'
#' @export
#' @examples
#' results <- load_example_results(multi_dataset=TRUE)
#' out <- plot_bar_dendro(results = results)
plot_bar_dendro <- function(results = load_example_results(multi_dataset = TRUE),
                            celltype_col = "cl_name",
                            target_branches = get_target_branches(),
                            keep_ancestors=names(target_branches),
                            facets = "ancestor_name",
                            add_test_target_celltypes=TRUE,
                            preferred_palettes = "tol",
                            legend.position="none",
                            heights = c(.3,1,.3,.3),
                            expand_dendro_x =rep(0.01,2),
                            q_threshold=0.05,
                            show_plot=TRUE,
                            save_path=NULL,
                            height = 16,
                            width = 13,
                            ...) {
  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  requireNamespace("ggdendro")
  requireNamespace("scales")
  # results=data.table::fread("/Users/bms20/Desktop/Rare Disease Celltyping/rare_disease_celltyping/results/phenomix_results.tsv.gz")
  hpo_id <- p <- fold_change <- ancestor_name <- ancestor_color <- cl_id <-
    enriched_phenotypes_norm <- enriched_phenotypes <- top_ancestor_name <-
    dummy <- NULL;

  {
    results <- results[,total_phenotypes:=data.table::uniqueN(hpo_id)]
    results <- HPOExplorer::add_hpo_name(results)
    results <- HPOExplorer::add_ancestor(results)
    results <- map_celltype(results)
    results[, enriched_phenotypes:=data.table::uniqueN(hpo_id[q<q_threshold],
                                                       na.rm = TRUE),
            by=c(celltype_col,"cl_id","ancestor","ancestor_name")]
    results[, phenotypes_per_ancestor:=data.table::uniqueN(hpo_id),
            by=c("ancestor","ancestor_name")]
    results_full <- data.table::copy(results)
    results <- KGExplorer::filter_dt(results,
                                     filters = list(ancestor_name=keep_ancestors))
    target_celltypes <- get_target_celltypes(target_branches=target_branches)
  }
  ## Convert ontology to dendrogram ontology
  {
    cl <- KGExplorer::get_ontology(name = "cl",
                                   add_ancestors = 1)
    ddata <- ontology_to_ggdendro(ont=cl,
                                  terms=as.character(unique(results$cl_id))
    )
    #### Make celltypes an ordered factor based on the clustering #####
    results[[celltype_col]] <- factor(results[[celltype_col]],
                                      levels=unique(ddata$labels$label),
                                      ordered = TRUE)
  }
  #### Filter the results ####
  by <- unique(c(celltype_col,"cl_id","ancestor","ancestor_name"))
  dat <- results[q<q_threshold & cl_id %in% unique(ddata$labels$id),
                 lapply(.SD,mean),
                 by=by,
                 .SDcols = c("enriched_phenotypes",
                             "phenotypes_per_ancestor",
                             "p","q","fold_change")
                 ]
  #### Get the top HPO category for each cell type ####
  add_top_value(dat=dat,
                sort_var="enriched_phenotypes",
                label_var=facets,
                group_var=celltype_col,
                new_var="top_ancestor_name",
                normalise_group=TRUE)
  #### Color each cell type x-axis label by the most commonly enriched HPO category ####
  cmap <- get_color_map(dat=dat,
                        columns = "top_ancestor_name",
                        ddata=ddata,
                        celltype_col=celltype_col,
                        preferred_palettes=preferred_palettes)
  color_map <- cmap$color_map
  color_vector <- cmap$color_vector

  #### Create summary bar plot ####
  #### Create tissue-celltype heatmap #####
  results[,dummy:=paste0("All"," (n=",total_phenotypes," phenotypes)")]
  tissue_plots <- plot_tissues(results = results,
                               facet_var = "dummy",
                               types = "bar")
  ggsummary <- tissue_plots$bar_plot +
                  ggplot2::labs(x=NULL, y=NULL) +
                  ggplot2::scale_fill_gradient(low="black",high="grey") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "transparent")
                   )
  #### Create faceted bar plot ####
  ggbars_out <- plot_bar_branches(results=dat,
                                  results_full=results_full,
                                  target_branches=target_branches,
                                  target_celltypes=target_celltypes,
                                  celltype_col=celltype_col,
                                  add_test_target_celltypes=add_test_target_celltypes,
                                  color_map=color_map,
                                  color_vector=color_vector,
                                  legend.position=legend.position,
                                  q_threshold=q_threshold,
                                  facets=facets)
  ggbars <- ggbars_out$plot
  #### Create dendrogram plot ####
  ggdend <-
    ggdendro::ggdendrogram(ddata) +
    ggplot2::scale_x_discrete(expand = expand_dendro_x) +
    ggdendro::theme_dendro() +
    ggplot2::scale_y_reverse()

  #### Plot -log(p) vs. N enrichments in a given target cell type
  ggprop_res <- plot_proportional_enrichment(results = results_full,
                                             color_map = color_map,
                                             target_branches=target_branches,
                                             show_plot = FALSE)

  #### Combine plots ####
  plot.margin <- ggplot2::unit(c(0,0,0,0),"cm")
  ggp <- patchwork::wrap_plots(ggsummary,
                               ggbars,
                               ggdend,
                               ggprop_res$plot,
                               ncol = 1,
                               heights = heights) +
    patchwork::plot_annotation(tag_levels = letters) &
    theme(plot.margin = plot.margin)
  #### Show plot ####
  if(show_plot) methods::show(ggp)
  #### Save plot ####
  if(!is.null(save_path)){
    KGExplorer::plot_save(plt = ggp,
                          path = save_path,
                          height = height,
                          width = width)
  }
  #### Return ####
  return(
    list(data=dat,
         plot=ggp)
  )
}

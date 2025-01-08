#' Plot bar dendrogram
#'
#' Create a plot summarising \pkg{MSTExplorer} results as a bar chart
#' with multiple facets showing selected branches from the
#' Human Phenotype Ontology (HPO). Also shows a dendrogram of celltype-celltype
#' relationships using the Cell Ontology (CL).
#' @param celltype_col Name of the cell type column in the \code{results}.
#' @param target_branches A named list of HPO
#' branches each matched with CL cell type branches that
#' correspond to on-target cell types across the two ontologies.
#' @param keep_ancestors Only HPO terms that have these ancestors will be kept.
#' @param add_test_target_celltypes Using the significant phenotype-cell type
#' association \code{results}, run proportional enrichment tests to
#' determine whether each cell type is overrepresented in a given HPO branch
#' relative to all other HPO branches. Overrepresented cell types will be
#' denoted by "*" above its bar.
#' @param expand_dendro_x Passed to \link[ggplot2]{scale_x_discrete}
#' in the cell type dendrogram.
#' @param cl Cell Ontology (CL) object from
#'  \code{KGExplorer::get_ontology("cl")}.
#' @inheritParams plot_
#' @inheritParams ggnetwork_plot_full
#' @inheritParams HPOExplorer::add_ancestor
#' @inheritParams KGExplorer::map_colors
#' @inheritParams KGExplorer::plot_save
#' @inheritParams ggplot2::facet_wrap
#' @inheritParams ggplot2::theme
#' @inheritDotParams EWCE::ewce_plot
#' @returns A bar chart with dendrogram of EWCE results in each cell type.
#'
#' @export
#' @examples
#' results <- load_example_results()
#' out <- plot_bar_dendro(results = results)
plot_bar_dendro <- function(results = load_example_results(),
                            celltype_col = "cl_name",
                            target_branches = get_target_branches(),
                            keep_ancestors=names(target_branches),
                            hpo = HPOExplorer::get_hpo(),
                            cl = KGExplorer::get_ontology(name = "cl",
                                                          lvl = 1,
                                                          remove_rings = TRUE),
                            facets = "ancestor_name",
                            add_test_target_celltypes=TRUE,
                            add_prop_test=FALSE,
                            preferred_palettes = "tol",
                            legend.position="none",
                            heights = c(.3,1,.15,.3),
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
  hpo_id <- cl_id <- sig_phenotypes <-
    dummy <- total_phenotypes <- phenotypes_per_ancestor <- NULL;

  {
    results <- results[,total_phenotypes:=data.table::uniqueN(hpo_id)]
    results <- HPOExplorer::add_hpo_name(results,
                                         hpo = hpo)
    results <- HPOExplorer::add_ancestor(results,
                                         hpo = hpo)
    results <- map_celltype(results)
    results[, sig_phenotypes:=data.table::uniqueN(hpo_id[q<q_threshold],
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
    ddata <- ontology_to_ggdendro(ont=cl,
                                  terms=as.character(unique(results$cl_id)))
    #### Make celltypes an ordered factor based on the clustering #####
    results[[celltype_col]] <- factor(results[[celltype_col]],
                                      levels=unique(ddata$labels$label),
                                      ordered = TRUE)
  }
  #### Filter the results ####
  by <- unique(c(celltype_col,"cl_id","ancestor","ancestor_name"))
  add_logfc(results)
  dat <- results[q<q_threshold & cl_id %in% unique(ddata$labels$id),
                 lapply(.SD,mean),
                 by=by,
                 .SDcols = c("sig_phenotypes",
                             "phenotypes_per_ancestor",
                             "p","q","logFC")
                 ]
  #### Get the top HPO category for each cell type ####
  add_top_value(dat=dat,
                sort_var="sig_phenotypes",
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
  ggbars_out <- plot_bar_dendro_facets(results=dat,
                                       results_full=results_full,
                                       target_branches=target_branches,
                                       target_celltypes=target_celltypes,
                                       celltype_col=celltype_col,
                                       add_test_target_celltypes=add_test_target_celltypes,
                                       color_map=color_map,
                                       color_vector=color_vector,
                                       legend.position=legend.position,
                                       q_threshold=q_threshold,
                                       add_prop_test=add_prop_test,
                                       facets=facets,
                                       hpo=hpo)
  ggbars <- ggbars_out$plot
  #### Create dendrogram plot ####
  ggdend <-
    ggdendro::ggdendrogram(ddata) +
    ggplot2::scale_x_discrete(expand = expand_dendro_x) +
    ggdendro::theme_dendro() +
    ggplot2::scale_y_reverse()

  #### Plot -log(p) vs. N enrichments in a given target cell type
  ggprop_out <- plot_proportional_enrichment(results = results_full,
                                             color_map = color_map,
                                             target_branches=target_branches,
                                             show_plot = FALSE)

  #### Combine plots ####
  plot.margin <- ggplot2::unit(c(0,0,0,0),"cm")
  ggp <- patchwork::wrap_plots(ggsummary,
                               ggbars,
                               ggdend,
                               ggprop_out$plot,
                               ncol = 1,
                               heights = heights) +
    patchwork::plot_annotation(tag_levels = letters) &
    ggplot2::theme(plot.margin = plot.margin)
  #### Show plot ####
  if(isTRUE(show_plot)) methods::show(ggp)
  #### Save plot ####
  if(!is.null(save_path)){
    KGExplorer::plot_save(plt = ggp,
                          save_path = save_path,
                          height = height,
                          width = width)
  }
  #### Return ####
  return(
    list(data=dat,
         plot=ggp,
         ggbars_out=ggbars_out,
         ggprop_out=ggprop_out)
  )
}

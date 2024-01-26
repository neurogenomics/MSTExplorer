#' EWCE plot
#'
#' A shallow wrapper for the \link[EWCE]{ewce_plot} function.
#' @inheritParams EWCE::ewce_plot
#' @inheritParams HPOExplorer::add_ancestor
#' @inheritParams ggnetwork_plot_full
#' @inheritDotParams EWCE::ewce_plot
#' @returns A bar chart with dendrogram of EWCE results in each cell type.
#'
#' @export
#' @importFrom EWCE ewce_plot
#' @examples
#' results <- load_example_results(multi_dataset=TRUE)
#' out <- multiewce_plot(results = results)
multiewce_plot <- function(results = load_example_results(multi_dataset = TRUE),
                           ancestor_names=c(
                             "Abnormality of the nervous system",
                             "Abnormality of the cardiovascular system",
                             "Abnormality of the immune system",
                             "Abnormality of the respiratory system",
                             "Abnormality of the eye",
                             "Abnormality of the endocrine system",
                             "Neoplasm"),
                      ...) {
  requireNamespace("ggplot2")

  results <- map_celltype(results)
  #### Get celltype dendrogram ####
  ## From CTD
  # ctd1 <- get_data("ctd_HumanCellLandscape.rds")
  # ctd2 <- get_data("ctd_DescartesHuman.rds")
  # dend_list <- EWCE:::prep_dendro(ctdIN = ctd[[1]])
  ## From ontology
  cl <- KGExplorer::get_ontology("cl", add_ancestors = 1)
  #### Method 1 ####
  # X <- KGExplorer::ontology_to(ont = cl,
  #                              to="similarity")
  # X <- X[unique(results$cell_type_ontology_term_id),
  #        unique(results$cell_type_ontology_term_id)]
  # hc <- stats::hclust(stats::dist(X) )
  # ddata <- ggdendro::dendro_data(hc)
  #### Method 2 ####
  # g <- KGExplorer::ontology_to(ont = cl,
  #                               to="tbl_graph")
  # KGExplorer::filter_graph(g, node_filters = list(name=results$cell_type_ontology_term_id))

  #### Method 3 ####
  hc <- KGExplorer::ontology_to(ont = cl,
                                to="igraph_dist_hclust")
  dend <- as.dendrogram(hc)
  dend2 <- dendextend::prune(dend,
                             leaves=labels(dend)[!labels(dend) %in% results$cell_type_ontology_term_id]
                             )
  ddata <- ggdendro::dendro_data(dend2)


  ddata$labels$label <- KGExplorer::map_ontology_terms(ont = cl,
                                                      terms=ddata$labels$label)
  ggdend <-
    ggdendro::ggdendrogram(ddata) +
    ggplot2::scale_x_discrete(expand =c(.01, .01)) +
    ggdendro::theme_dendro() +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0,0,0,0), "cm"))
  # ggplot(ggdendro::segment(ddata)) +
  #   ggplot2::geom_segment(ggplot2::aes_string(
  #     x = "x", y = "y",
  #     xend = "xend", yend = "yend"
  #   )) +
  #   ggdendro::theme_dendro()
  # b1 <- ggplot(ggdendro::segment(ddata)) +
  #   ggplot2::geom_segment(ggplot2::aes(
  #     x = x, y = y,
  #     xend = xend, yend = yend
  #   )) +
  #   ggdendro::theme_dendro()
  # EWCE::ewce_plot(total_res = results,
  #                 mtc_method = mtc_method,
  #                 ctd = ctd,
  #                 ...)
  results <- HPOExplorer::add_ancestor(phenos = results)
  dat <- results[q<0.05, list(enriched_phenotypes=data.table::uniqueN(hpo_id),
                              p=mean(p),
                              q=mean(q),
                              fold_change=mean(fold_change)),
                 by=c("CellType","cell_type","cell_type_ontology_term_id",
                      "ancestor","ancestor_name")]
  dat <- dat[cell_type %in% unique(ddata$labels$label),]
  dat <- dat[ancestor_name %in% ancestor_names,]

  dat$cell_type <- factor(dat$cell_type,
                          levels=unique(ddata$labels$label),
                          ordered = TRUE)
  #### Get the top HPO category for each cell type ####
  ### NOrmalise count within each ancestor
  dat[,enriched_phenotypes_norm:=enriched_phenotypes/max(enriched_phenotypes),
      by=c("ancestor")]
  dat[,top_ancestor_name:=head(ancestor_name[which.max(enriched_phenotypes_norm)],1),
      by=c("cell_type","cell_type_ontology_term_id")]
  # cellmeta <- data.table::data.table(
  #   cl@elementMetadata[,c("name","ancestor","ancestor_name")],
  #   key = "name"
  # )[unique(ddata$labels$label),]

  color_dict <- KGExplorer::map_colors(dat,
                                         columns = "top_ancestor_name",
                                         preferred_palettes = "brewer.set2",
                                         as="dict")[[1]]
  dat2 <- unique(dat[,c("cell_type","top_ancestor_name")]
                 )[,color:=color_dict[top_ancestor_name]] |>
    data.table::setkeyv("cell_type")
  color_vector <- dat2[unique(ddata$labels$label),]$color
  color_vector[is.na(color_vector)] <- "grey20"
  # color_vector <- KGExplorer::map_colors(unique(dat[,c("cell_type","top_ancestor_name")]),
  #                                        columns = "top_ancestor_name",
  #                                        preferred_palettes = "brewer.set2",
  #                                      as="vector")[[1]]

  ggbars <- ggplot2::ggplot(dat, ggplot2::aes(x=cell_type,
                                              y=enriched_phenotypes,
                                              fill=ancestor_name)) +
    ggplot2::geom_bar(stat="identity", show.legend = FALSE) +
    ggplot2::scale_fill_manual(values=color_dict) +
    ggplot2::facet_wrap(facets = "ancestor_name",
                        scales="free_y", ncol = 1) +
    ggplot2::labs(x=NULL,y="Enriched phenotypes") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90, hjust = 1, vjust=0.5,
      color = unname(color_vector)
      ),
          strip.background = ggplot2::element_rect(fill = "transparent")

          )

  patchwork::wrap_plots(ggbars ,
                        ggdend + ggplot2::scale_y_reverse(),

                        ncol = 1, heights = c(1,.3))

}

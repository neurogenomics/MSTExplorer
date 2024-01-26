#' Plot bar dendrogram
#'
#' Create a plot summarising MultiEWCE results as a bar chart with multiple
#'  facets and a cell ontology-based dendrogram.
#' @param ancestor_names Keep terms with specific ancestors.
#' @param celltype_col Name of the cell type column in the \code{results}.
#' @inheritParams ggnetwork_plot_full
#' @inheritParams HPOExplorer::add_ancestor
#' @inheritParams ggplot2::facet_wrap
#' @inheritDotParams EWCE::ewce_plot
#' @returns A bar chart with dendrogram of EWCE results in each cell type.
#'
#' @export
#' @examples
#' results <- load_example_results(multi_dataset=TRUE)
#' out <- plot_bar_dendro(results = results)
plot_bar_dendro <- function(results = load_example_results(multi_dataset = TRUE),
                            celltype_col = "cl_name",
                            ancestor_names=c(
                              "Abnormality of the nervous system",
                              "Abnormality of the cardiovascular system",
                              "Abnormality of the immune system",
                              "Abnormality of the respiratory system",
                              "Abnormality of the eye",
                              "Abnormality of the endocrine system",
                              "Neoplasm"),
                            facets = "ancestor_name",
                            ...) {
  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  requireNamespace("dendextend")
  requireNamespace("ggdendro")
  hpo_id <- p <- fold_change <- cl_name <- ancestor_name <-
    enriched_phenotypes_norm <- enriched_phenotypes <- top_ancestor_name <-
    color <- NULL;

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
  # X <- X[unique(results$cl_id),
  #        unique(results$cl_id)]
  # hc <- stats::hclust(stats::dist(X) )
  # ddata <- ggdendro::dendro_data(hc)
  #### Method 2 ####
  # g <- KGExplorer::ontology_to(ont = cl,
  #                               to="tbl_graph")
  # KGExplorer::filter_graph(g,
  # node_filters = list(name=results$cl_id))

  #### Method 3 ####
  hc <- KGExplorer::ontology_to(ont = cl,
                                to="igraph_dist_hclust")
  dend <- stats::as.dendrogram(hc)
  messager("Pruning dendrogram.")
  dend2 <- dendextend::prune(
    dend,
    leaves=labels(dend)[!labels(dend) %in% results$cl_id])
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
                 by=c(celltype_col,"ancestor","ancestor_name")]

  dat <- dat[get(celltype_col) %in% unique(ddata$labels$label),]
  dat <- KGExplorer::filter_dt(dat,
                               filters = list(ancestor_name=ancestor_names))

  dat[[celltype_col]] <- factor(dat[[celltype_col]],
                                levels=unique(ddata$labels$label),
                                ordered = TRUE)
  #### Get the top HPO category for each cell type ####
  ### NOrmalise count within each ancestor
  dat[,enriched_phenotypes_norm:=(enriched_phenotypes/max(enriched_phenotypes)),
      by=c("ancestor")]
  dat[,top_ancestor_name:=head(ancestor_name[which.max(enriched_phenotypes_norm)],1),
      by=c(celltype_col)]
  # cellmeta <- data.table::data.table(
  #   cl@elementMetadata[,c("name","ancestor","ancestor_name")],
  #   key = "name"
  # )[unique(ddata$labels$label),]

  color_dict <- KGExplorer::map_colors(dat,
                                         columns = "top_ancestor_name",
                                         preferred_palettes = "brewer.set2",
                                         as="dict")[[1]]
  dat2 <- unique(dat[,c(celltype_col,"top_ancestor_name"), with=FALSE]
                 )[,color:=color_dict[top_ancestor_name]] |>
    data.table::setkeyv(celltype_col)
  color_vector <- dat2[unique(ddata$labels$label),]$color
  color_vector[is.na(color_vector)] <- "grey20"
  # color_vector <- KGExplorer::map_colors(unique(dat[,c(celltype_col,"top_ancestor_name")]),
  #                                        columns = "top_ancestor_name",
  #                                        preferred_palettes = "brewer.set2",
  #                                      as="vector")[[1]]

  ggbars <- ggplot2::ggplot(dat, ggplot2::aes(x=!!ggplot2::sym(celltype_col),
                                              y=enriched_phenotypes,
                                              fill=ancestor_name)) +
    ggplot2::geom_bar(stat="identity", show.legend = FALSE) +
    ggplot2::scale_fill_manual(values=color_dict) +
    ggplot2::facet_wrap(facets = facets,
                        scales="free_y", ncol = 1) +
    ggplot2::labs(x=NULL,y="Enriched phenotypes") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90, hjust = 1, vjust=0.5,
      color = unname(color_vector)
      ),
      strip.background = ggplot2::element_rect(fill = "transparent")
      )

  ggp <- patchwork::wrap_plots(ggbars ,
                        ggdend + ggplot2::scale_y_reverse(),
                        ncol = 1, heights = c(1,.3))
  #### Return ####
  return(
    list(data=dat,
         plot=ggp)
  )
}

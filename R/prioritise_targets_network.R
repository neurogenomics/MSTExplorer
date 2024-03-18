#' Prioritise targets network
#'
#' Plot the output of \link[MSTExplorer]{prioritise_targets}
#' as an interactive network with \link[visNetwork]{visNetwork}.
#' See \link[igraph]{layout} for a list of all layout functions
#' that can be passed to the \code{prioritise_targets_network(layout=)} argument
#' as a string (e.g. "layout.fruchterman.reingold" or "layout_as_tree").
#' See \href{https://rpubs.com/JanpuHou/337696}{here} for some visualized
#'  examples of each layout.
#' @source \href{https://book.archnetworks.net/visualization}{
#' book.archnetworks.net}
#' @source
#' \href{https://r-graph-gallery.com/310-custom-hierarchical-edge-bundling}{
#' Edge bundling with \pkg{ggraph}
#' }
#' @source \href{https://visjs.github.io/vis-network/examples/}{
#' visNetwork examples}
#' @param top_targets output of \link[MSTExplorer]{prioritise_targets}.
#' @param vertex_vars Columns within \code{top_targets}
#' to include as vertices/nodes within the network.
#' @param group_var Variable to group nodes by.
#' @param edge_color_var Variable to color edges by.
#' @param edge_size_var Variable to scale edges by.
#' @param mediator_var Variable to connect cell types and phenotypes by
#'  (i.e. "gene_symbol").
#'  If \code{NULL}, instead will connect node hierarchically:
#'  ancestor_name --> Phenotype --> CellType --> gene_symbol
#' @param degree of depth of nodes to be colored. Default to 1.
#' Set high number to have the entire sub-network.
#'  In case of "hierarchical" algorithm, you can also pass a
#'  list(from = 1, to = 1) to control degree in both direction.
#' @param verbose Print messages.
#' @inheritParams prioritise_targets
#' @inheritParams KGExplorer::plot_graph_visnetwork
#' @returns A named list containing the \link[visNetwork]{visNetwork} plot
#' and the the graph used to make the plot.
#'
#' @export
#' @examples
#' top_targets <- MSTExplorer::example_targets$top_targets[1:10]
#' top_targets[,estimate:=fold_change]
#' top_targets <- map_celltype(top_targets)
#' vn <- prioritise_targets_network(top_targets = top_targets)
prioritise_targets_network <- function(top_targets,
                                       vertex_vars = c("disease_name",
                                                       # "ancestor_name",
                                                       "hpo_name",
                                                       "cl_name",
                                                       "gene_symbol"),
                                       group_var = vertex_vars[[1]],
                                       colour_var = "node_type",
                                       edge_color_var = "logFC",
                                       edge_size_var = "logFC",
                                       # mediator_var = list(),
                                       mediator_var = list(c(1,2),c(2,3),c(3,4),c(1,3)),
                                       layout = "layout_with_sugiyama",
                                       solver = "forceAtlas2Based",
                                       physics = FALSE,
                                       forceAtlas2Based = list(
                                         avoidOverlap=.5,
                                         gravitationalConstant=-50),
                                       scaling=NULL,
                                       arrows = "from",
                                       smooth=list(enabled=TRUE,
                                                   type="cubicBezier",
                                                   roundness=.5),
                                       add_visExport = FALSE,
                                       degree = 1,
                                       width = "100%",
                                       height = "90vh",
                                       main = "Rare Disease Celltyping",
                                       submain = "Prioritised Targets Network",
                                       preferred_palettes = "kovesi.linear_bmy_10_95_c78",
                                       randomSeed = 2023,
                                       verbose = TRUE,
                                       show_plot = TRUE,
                                       save_path = tempfile(
                                         fileext =
                                           "_prioritise_targets_network.html"
                                       ),
                                       run_prune_ancestors=FALSE
                                       ){
  requireNamespace("ggplot2")
  requireNamespace("pals")

  add_logfc(top_targets)
  if(any(c("cl_name","cl_id") %in% vertex_vars)){
    map_celltype(top_targets)
  }
  if(isTRUE(run_prune_ancestors)){
    hpo <- HPOExplorer::get_hpo()
    top_targets <- KGExplorer::prune_ancestors(dat = top_targets,
                                               id_col = "hpo_id",
                                               ont = hpo)
  }
  #### Network ####
  g <- targets_to_graph(top_targets = top_targets,
                        vertex_vars = c(group_var,vertex_vars),
                        group_var = group_var,
                        edge_color_var = edge_color_var,
                        edge_size_var = edge_size_var,
                        mediator_var = mediator_var,
                        verbose = verbose)
  out <- KGExplorer::plot_graph_visnetwork(
    g = g,
    label_var = "name",
    colour_var = colour_var,
    save_path = save_path,
    preferred_palettes = preferred_palettes,
    layout = layout,
    solver = solver,
    physics = physics,
    forceAtlas2Based = forceAtlas2Based,
    scaling = scaling,
    arrows = arrows,
    smooth = smooth,
    add_visExport = add_visExport,
    degree = degree,
    height = height,
    width = width,
    main = main,
    submain = submain,
    randomSeed = randomSeed,
    show_plot = show_plot)
  return(out)
}



# edgebundle::install_bundle_py()
# reticulate::use_condaenv(condaenv = "r-reticulate")
# hbundle <-   edgebundle::edge_bundle_hammer(object = g,
#                                             xy =  igraph::layout_with_kk(g),
#                                             bw = 5,
#                                             decay = 0.3)
# #
# f6_3c <-   ggplot() +
#   geom_path(data = hbundle, aes(x, y, group = group),
#             col = "gray66", size = 0.5) +
#   geom_point(data = xy, aes(x, y, col = Region),
#              size = 5, alpha = 0.75, show.legend = FALSE) +
#   theme_void()

#   ggraph::ggraph(g,
#                  layout = 'dendrogram') +
#     ggraph::geom_edge_density() +
#     ggraph::geom_edge_bend() +
#     ggraph::geom_node_point(ggplot2::aes(color=group, shape=shape))
#     # ggraph::geom_node_label(ggplot2::aes(label=name, fill=group))
#     ggraph::geom_conn_bundle(data = ggraph::get_con(from = edges$from,
#                                                     to = edges$to),
#                              alpha=1,
#                              colour="black",
#                              tension = .5) +
#     ggraph::geom_node_point() +
#     ggplot2::theme_void()



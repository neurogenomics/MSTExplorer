#' Prioritise targets network
#'
#' Plot the output of \link[MultiEWCE]{prioritise_targets}
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
#' @param top_targets output of \link[MultiEWCE]{prioritise_targets}.
#' @param vertex_vars Columns within \code{top_targets}
#' to include as vertices/nodes within the network.
#' @param group_var Variable to group nodes by.
#' @param edge_color_var Variable to color edges by.
#' @param edge_size_var Variable to scale edges by.
#' @param mediator_var Variable to connect cell types and phenotypes by
#'  (i.e. "Gene").
#'  If \code{NULL}, instead will connect node hierarchically:
#'  ancestor_name --> Phenotype --> CellType --> Gene
#' @param degree of depth of nodes to be colored. Default to 1.
#' Set high number to have the entire sub-network.
#'  In case of "hierarchical" algorithm, you can also pass a
#'  list(from = 1, to = 1) to control degree in both direction.
#' @param show_plot Print the plot after it has been created.
#' @param save_path Path to save HTML version of network to.
#' @param verbose Print messages.
#' @param add_visExport Add PDF download button.
#'
#' @inheritParams visNetwork::visIgraph
#' @inheritParams visNetwork::visIgraphLayout
#' @inheritParams visNetwork::visPhysics
#' @inheritParams visNetwork::visNodes
#' @inheritParams visNetwork::visEdges
#' @inheritParams visNetwork::visOptions
#' @inheritParams visNetwork::visNetwork
#' @returns A named list containing the \link[visNetwork]{visNetwork} plot
#' and the the graph used to make the plot.
#'
#' @export
#' @import ggplot2
#' @importFrom stats setNames
#' @importFrom dplyr %>%
#' @examples
#' top_targets <- MultiEWCE::example_targets$top_targets
#' vn <- prioritise_targets_network(top_targets = top_targets)
prioritise_targets_network <- function(top_targets,
                                       vertex_vars = c("DiseaseName",
                                                       "ancestor_name",
                                                       "Phenotype",
                                                       "CellType",
                                                       "Gene"),
                                       group_var = vertex_vars[[1]],
                                       edge_color_var = group_var,
                                       edge_size_var = "fold_change",
                                       mediator_var = list(),
                                       show_plot = TRUE,
                                       save_path = tempfile(
                                         fileext =
                                           "_prioritise_targets_network.html"
                                         ),
                                       layout = "layout_with_kk",
                                       solver = "forceAtlas2Based",
                                       physics = FALSE,
                                       forceAtlas2Based = list(
                                         avoidOverlap=.5,
                                         gravitationalConstant=-50),
                                       scaling=list(label=
                                                      list(drawThreshold=0)
                                                    ),
                                       smooth=list(enabled=TRUE,
                                                   type="cubicBezier",
                                                   roundness=.5),
                                       add_visExport = FALSE,
                                       degree = if(is.null(mediator_var)){
                                         1
                                       } else if (is.list(mediator_var)){
                                         if(length(mediator_var)==0){
                                           2
                                         } else {
                                           1
                                         }
                                       },
                                       height = NULL,
                                       width = NULL,
                                       main = "Rare Disease Celltyping",
                                       submain = "Prioritised Targets Network",
                                       randomSeed = 2023,
                                       verbose = TRUE
                                       ){
  # devoptera::args2vars(prioritise_targets_network)

  requireNamespace("ggplot2")
  requireNamespace("pals")

  #### Network ####
  g <- targets_to_graph(top_targets = top_targets,
                        vertex_vars = c(group_var,vertex_vars),
                        group_var = group_var,
                        edge_color_var = edge_color_var,
                        edge_size_var = edge_size_var,
                        mediator_var = mediator_var,
                        node_palette = pals::kovesi.linear_bmy_10_95_c78,
                        # format = "ggnetwork",
                        verbose = verbose)

  # plt <- HPOExplorer::network_3d(g = g,
  #                                node_color_var = "node_type",
  #                                node_symbol_var = "node_type",
  #                                add_labels = FALSE)

  plt <- plot_visnetwork(g = g,
                         save_path = save_path,
                         layout = layout,
                         solver = solver,
                         physics = physics,
                         forceAtlas2Based = forceAtlas2Based,
                         scaling = scaling,
                         smooth = smooth,
                         add_visExport = add_visExport,
                         degree = degree,
                         height = height,
                         width = width,
                         main = main,
                         submain = submain,
                         randomSeed = randomSeed,
                         verbose = verbose)
  #### Show ####
  if(isTRUE(show_plot)) {
    # utils::browseURL(save_path)
    methods::show(plt)
  }
  return(list(plot=plt,
              graph=g))


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



}

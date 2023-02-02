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
#' @param top_targets output of \link[MultiEWCE]{prioritise_targets}.
#' @param vertex_vars Columns within \code{top_targets}
#' to include as vertices/nodes within the network.
#' @param group_var Variable to group nodes by.
#' @param edge_var Variable to color edges by.
#' @param show_plot Print the plot after it has been created.
#' @param save_path Path to save HTML version of network to.
#' @param verbose Print messages.
#'
#' @inheritParams visNetwork::visIgraph
#' @inheritParams visNetwork::visIgraphLayout
#' @returns A named list containing the \link[visNetwork]{visNetwork} plot
#' and the the graph used to make the plot.
#'
#' @export
#' @import ggplot2
#' @importFrom stats setNames
#' @examples
#' top_targets <- prioritise_targets()
#' vn <- prioritise_targets_network(top_targets = top_targets)
prioritise_targets_network <- function(top_targets,
                                       vertex_vars = c("Phenotype",
                                                       "CellType",
                                                       "Gene"),
                                       group_var = "ancestor_name",
                                       edge_var = "fold_change",
                                       show_plot = TRUE,
                                       save_path = tempfile(
                                         fileext =
                                           "_prioritise_targets_network.html"
                                         ),
                                       layout = "layout_with_kk",
                                       randomSeed = 2023,
                                       verbose = TRUE
                                       ){
  # templateR:::source_all()
  # templateR:::args2vars(prioritise_targets_network)

  requireNamespace("ggplot2")
  requireNamespace("visNetwork")

  #### Network ####
  g <- targets_to_graph(top_targets = top_targets,
                        vertex_vars = c(group_var,vertex_vars),
                        group_var = group_var,
                        edge_var = edge_var,
                        verbose = verbose)
  #### Plot ####
  messager("Creating plot.",v=verbose)
  visnet <- visNetwork::visIgraph(g, randomSeed = randomSeed) |>
    visNetwork::visIgraphLayout(layout = layout) |>
    visNetwork::visNodes(font = list(color="#F0FFFF",
                                     strokeWidth=2,
                                     strokeColor="rgba(0,0,0,0.8)"
                                     ),
                         shadow = list("enabled"=TRUE),
                         opacity = 0.9) |>
    visNetwork::visEdges(shadow = list("enabled"=FALSE),
                         # smooth = list("enabled"=TRUE,
                         #               type="straightCross",
                         #               roundness=.5),
                         color = list(opacity = 0.9)) |>
    visNetwork::visLegend() |>
    # visNetwork::visGroups(groupname = unique(igraph::vertex_attr(g,"group"))[[1]],
    #                       color="brown")
    # visNetwork::visClusteringByGroup(groups = igraph::vertex_attr(g,"group"))
    visNetwork::visExport(type = "pdf") |>
    visNetwork::visOptions(height=700, width = 1300)
  if(!is.null(save_path)) {
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    messager("Saving plot ==>",save_path,v=verbose)
    visNetwork::visSave(visnet,file = save_path)
    # {
    #   grDevices::pdf(file = "~/Downloads/network2.pdf")
    #   methods::show(visnet)
    #   grDevices::dev.off()
    # }
  }
  #### Show ####
  if(isTRUE(show_plot)) {
    # utils::browseURL(save_path)
    methods::show(visnet)
  }
  return(list(plot=visnet,
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


  # n <- intergraph::asNetwork(g)
  # ggplot(n,
  #                        aes(x = x,
  #                            y = y,
  #                            xend = xend,
  #                            yend = yend)) +
  #   # geom_label(aes_string(label = "vertex.names", colour = "group")) +
  #   geom_point(aes(color=node_type, shape=node_type)) +
  #   geom_text(aes(label = node), color = "black") +
  #   ggnetwork::geom_edges(aes(color=fold_change)) +
  #   # scale_colour_gradient2(low = "white", mid = "yellow", high = "red") +
  #   scale_size(trans = "exp") +
  #   # guides(size = "none") +
  #   # labs(colour = colour_label) +
  #   theme_void()

#
#
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

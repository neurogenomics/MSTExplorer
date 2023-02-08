plot_visnetwork <- function(g,
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
                            height = 1000,
                            width = 1300,
                            randomSeed = 2023,
                            verbose=TRUE){
  requireNamespace("visNetwork")

  messager("Creating plot.",v=verbose)
  visnet <- visNetwork::visIgraph(g,
                                  randomSeed = randomSeed) |>
    visNetwork::visIgraphLayout(layout = layout,
                                physics = physics) |>
    visNetwork::visPhysics(solver=solver,
                           forceAtlas2Based=forceAtlas2Based,
                           enabled = physics) |>
    visNetwork::visNodes(font = list(color="#F0FFFF",
                                     strokeWidth=2,
                                     strokeColor="rgba(0,0,0,1)"
    ),
    shadow = list(enabled=TRUE,
                  size = 10),
    opacity = 0.75,
    color = list(hover=list(background="rgba(0,0,0,.5)")
    ),
    scaling = scaling
    ) |>
    visNetwork::visEdges(shadow = list(enabled=FALSE),
                         smooth = smooth,
                         color = list(opacity = 0.75)) |>
    # visNetwork::visLegend() |>
    # visNetwork::visClusteringByConnection(nodes = unique(top_targets[[group_var]])) |>
    # visNetwork::visGroups(groupname = unique(igraph::vertex_attr(g,"group"))[[2]],
    #                       color="green")
    # visNetwork::visClusteringByGroup(groups = igraph::vertex_attr(g,"group"))
    visNetwork::visInteraction(hover = TRUE) |>
    visNetwork::visOptions(height=height,
                           width = width,
                           highlightNearest = list(enabled=TRUE,
                                                   degree=2))
  if(add_visExport){
    visnet <-  visnet |> visNetwork::visExport(type = "pdf")
  }
  #### Save network ####
  if(!is.null(save_path)) {
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    messager("Saving plot ==>",save_path,v=verbose)
    visNetwork::visSave(visnet,
                        file = save_path,
                        selfcontained = TRUE,
                        background = "transparent")
    return(visnet)
  }
  #### Plot ####

  # {
  #   grDevices::pdf(file = "~/Downloads/network2.pdf")
  #   methods::show(visnet)
  #   grDevices::dev.off()
  # }
}

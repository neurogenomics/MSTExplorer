targets_to_graph <- function(top_targets,
                             vertex_vars,
                             group_var,
                             edge_color_var=group_var,
                             edge_size_var="fold_change",
                             mediator_var = "Gene",
                             format="visnetwork",
                             verbose=TRUE){

  requireNamespace("igraph")
  fold_change <- node_type <- shape <- color <-
    Phenotype <- ancestor_name <- NULL;

  messager("Creating network.",v=verbose)
  #### Create vertices ####
  vertex_vars <- unique(vertex_vars)
  shapes <- c("database","circle","box")
  if("ancestor_name" %in% vertex_vars && length(vertex_vars)==4){
    shapes <- c("database",shapes)
  }
  shape_dict <- stats::setNames(
    shapes,
    unique(vertex_vars))
  color_dict <- stats::setNames(
    pals::ocean.thermal(length(unique(vertex_vars))+2)[-1],
    unique(vertex_vars))
  ##### Remove Phenotypes that are also ancestor #####
  if("ancestor_name" %in% names(top_targets)){
    top_targets <- top_targets[Phenotype!=ancestor_name,]
  }
  #### Make vertex metadata ####
  vertices <- (
    data.table::melt.data.table(
      top_targets[,c(vertex_vars,group_var),with=FALSE],
      id.vars = group_var,
      measure.vars = vertex_vars,
      variable.name = "node_type",
      value.name = "node") |>
      rev() |> unique()
  )[,shape:=shape_dict[node_type]][,color:=color_dict[node_type]] |>
    data.table::setkeyv(cols = "node")
  #### ancestor_name is only relevant metadata for Phenotype nodes ####
  vertices[node_type!="Phenotype",]$ancestor_name <- NA
  vertices <- unique(vertices)
  #### Merge graphs ####
  if(!is.null(mediator_var)){
    subgraphs <- lapply(vertex_vars[vertex_vars!=mediator_var], function(v){
      vv <- c(v,mediator_var)
      dt <- unique(
        top_targets[,vv,with=FALSE]
      )
      igraph::graph_from_data_frame(dt)
    })
  } else{
    subgraphs <- lapply(seq_len(length(vertex_vars)-1), function(i){
      vv <- vertex_vars[c(i,i+1)]
      dt <- unique(
        top_targets[,c(vv,c(edge_color_var,edge_size_var)),
                    with=FALSE][,fold_change:=mean(fold_change), by=vv]
      )
      igraph::graph_from_data_frame(dt)
    })
  }

  g <- igraph::graph.union(subgraphs)
  #### Name edges ####
  # igraph::edge_attr(g,"id") <- paste0("edge",seq_len(length(igraph::E(g))))
  #### Add nodes metadata ####
  for(x in names(vertices)){
    igraph::vertex_attr(g,x) <- as.character(
      vertices[names(igraph::V(g))][[x]]
    )
  }
  if(group_var %in% names(vertices)){
    igraph::vertex_attr(g,"group") <- vertices[names(igraph::V(g))][[group_var]]
  }
  #### Add hoverdata ####
  igraph::vertex_attr(g,"title") <- lapply(seq_len(length(g)), function(i){
    nms <- igraph::vertex_attr_names(g)
    nms <- nms[!nms %in% c('shape','color','value','name')]
    lapply(nms, function(nm){
      value <- igraph::vertex_attr(g,nm)[i]
      if(!is.na(value)) {
        paste0("<strong>",nm,"</strong>: ",value)
      } else {
        ""
      }
    }) |>
      paste(collapse = "<br>") |>
      gsub(pattern="<br><br>",replacement="<br>")
  }) |> unlist()
  # lapply(stats::setNames(unique(vertices[[group_var]]),
  #                        unique(vertices[[group_var]])),
  #        function(group){
  #          seed_nodes <- vertices[get(group_var)==group,]$node
  #          # g2 <- igraph::induced_subgraph(graph = g, vids = group_nodes)
  #          adj_nodes <- igraph::adjacent_vertices(graph = g,
  #                                                 v = seed_nodes,
  #                                                 mode = "all")
  #          igraph::edge_attr(g)$id
  #          group_nodes <- unique(c(seed_nodes,
  #                                  unname(unlist(lapply(adj_nodes, names)))))
  #
  #        })
  igraph::vertex_attr(g,"name") <- stringr::str_wrap(
    vertices[names(igraph::V(g))][["node"]],
    width = 10)
  #### Add edge color ####
  # edge_color <- igraph::edge_attr(g)[grep(edge_color_var,igraph::edge_attr_names(g), value = TRUE)] |>
  #   data.table::as.data.table() |>
  #   data.table::fcoalesce()
  # edge_color_dict <- stats::setNames(
  #   pals::alphabet(length(unique(vertices[[group_var]]))),
  #   unique(vertices[[group_var]]))
  # igraph::edge_attr(g,"color") <- edge_color_dict[edge_color]
  #### Add edge size ####
  # edge_size <- igraph::edge_attr(g)[grep(edge_size_var,igraph::edge_attr_names(g), value = TRUE)] |>
  #   data.table::as.data.table() |>
  #   data.table::fcoalesce()
  # cols <- c(edge_size_var)
  # top_targets[,c(vertex_vars,edge_color_var,edge_size_var), with=FALSE][,(cols):=lapply(.SD,mean),.SDcols=cols,by=c("Phenotype")]

  # igraph::edge_attr(g,"width") <- edge_size
  # g <- igraph::simplify(g,)
  if(isTRUE(format=="ggnetwork")){
    requireNamespace("ggnetwork")
    g2 <- ggnetwork::fortify(g)
    rownames(g2) <- paste0("edge",seq_len(nrow(g2)))
    return(g2)
  } else{
    return(g)
  }
}

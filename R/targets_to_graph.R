targets_to_graph <- function(top_targets,
                             vertex_vars,
                             group_var,
                             edge_var,
                             verbose=TRUE){

  requireNamespace("igraph")
  fold_change <- node_type <- shape <- color <-
    Phenotype <- ancestor_name <- NULL;

  vertex_vars <- unique(vertex_vars)
  messager("Creating network.",v=verbose)
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
      top_targets[,vertex_vars,with=FALSE],
      measure.vars = vertex_vars,
      variable.name = "node_type",
      value.name = "node") |>
      rev() |> unique()
  )[,shape:=shape_dict[node_type]][,color:=color_dict[node_type]] |>
    data.table::setkeyv(cols = "node")

  #### Merge graphs ####
  subgraphs <- lapply(seq_len(length(vertex_vars)-1), function(i){
    vv <- vertex_vars[c(i,i+1)]
    dt <- unique(
      top_targets[,c(vv,edge_var), with=FALSE][,fold_change:=mean(fold_change),
                                               by=vv]
    )
    igraph::graph_from_data_frame(dt)
  })
  g <- do.call(igraph::union, subgraphs)
  #### Add nodes metadata ####
  for(x in names(vertices)){
    igraph::vertex_attr(g,x) <- as.character(
      vertices[names(igraph::V(g))][[x]]
    )
  }
  if(group_var %in% names(vertices)){
    igraph::vertex_attr(g,"group") <- vertices[names(igraph::V(g))][[group_var]]
  }
  igraph::vertex_attr(g,"name") <- stringr::str_wrap(
    vertices[names(igraph::V(g))][["node"]],
    width = 10)
  #### Add edge metadata ####
  edge_data <- igraph::edge_attr(g) |>
    data.table::as.data.table() |>
    data.table::fcoalesce()
  for(e in grep(edge_var,igraph::edge_attr_names(g), value = TRUE)){
    g <- igraph::remove.edge.attribute(g, e)
  }
  for(e in c(edge_var,"width","color")){
    igraph::edge_attr(g,e) <- edge_data
  }
  # g <- igraph::simplify(g,)

  # igraph::vertex.attributes(g)
  # igraph::edge.attributes(g)
  # visNetwork::visNetwork(g)
  return(g)
}

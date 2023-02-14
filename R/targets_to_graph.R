targets_to_graph <- function(top_targets,
                             vertex_vars,
                             group_var,
                             metadata_vars=c("HPO_ID","definition",
                                             "ontLvl","tier_merge"),
                             edge_color_var=group_var,
                             edge_size_var="fold_change",
                             mediator_var = "Gene",
                             node_palette = pals::isol, #pals::ocean.thermal,
                             edge_palette = node_palette,
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
  ##### Remove Phenotypes that are also ancestor #####
  if("ancestor_name" %in% names(top_targets)){
    top_targets <- top_targets[Phenotype!=ancestor_name,]
  }
  #### Make vertex metadata ####
  vertices <- (
    data.table::melt.data.table(
      top_targets[,c(vertex_vars,group_var,metadata_vars),with=FALSE],
      id.vars = c(group_var,metadata_vars),
      measure.vars = vertex_vars,
      variable.name = "node_type",
      value.name = "node") |>
      rev() |> unique()
  )
  #### Add node shapes ####
  shape_dict <- stats::setNames(
    shapes,
    unique(vertex_vars))
  vertices[,shape:=shape_dict[node_type]]
  #### Add node colors ####
  if(!is.null(node_palette)){
    color_dict <- stats::setNames(
      node_palette(length(unique(vertex_vars))+2)[-1],
      unique(vertex_vars))
    vertices[,color:=color_dict[node_type]] |>
      data.table::setkeyv(cols = "node")
  }
  vertices <- vertices[,utils::head(.SD, 1),by = c("node")]
  #### ancestor_name is only relevant metadata for Phenotype nodes ####
  vertices[node_type!="Phenotype",]$ancestor_name <- NA
  vertices <- unique(vertices)
  vertices$name <- stringr::str_wrap(vertices$node,
                                     width = 10)
  #### Merge graphs ####
  if(is.character(mediator_var)){
    subgraphs <- lapply(vertex_vars[vertex_vars!=mediator_var], function(v){
      vv <- c(v,mediator_var)
      dt <- unique(
        top_targets[,vv,with=FALSE]
      )
      igraph::graph_from_data_frame(dt)
    })
  } else if(is.list(mediator_var)){
    ilist <- if(length(mediator_var)==0){
      list(c(1,2),c(2,3),c(3,4),c(2,4))
    } else {
      mediator_var
    }
    if(length(unique(unlist(ilist)))>length(vertex_vars)){
      stp <- paste("When mediator_var is a list, mediator_var must be",
                   "equal to or less than the length of vertex_vars.")
      stop(stp)
    }
    subgraphs <- lapply(ilist, function(il){
      vv <- vertex_vars[il]
      cols <- edge_size_var
      dt <- unique(
        top_targets[,c(vv,c(edge_color_var,edge_size_var)),
                    with=FALSE][,(cols):=lapply(.SD,mean),.SDcols=cols,by=vv]
      )
      igraph::graph_from_data_frame(dt)
    })
  } else {
    subgraphs <- lapply(seq_len(length(vertex_vars)-1), function(i){
      vv <- vertex_vars[c(i,i+1)]
      dt <- unique(
        top_targets[,c(vv,c(edge_color_var,edge_size_var)),
                    with=FALSE][,fold_change:=mean(fold_change), by=vv]
      )
      igraph::graph_from_data_frame(dt)
    })
  }
  #### Merge subgraphs ####
  g <- igraph::graph.union(subgraphs)
  #### Name edges ####
  # igraph::edge_attr(g,"id") <- paste0("edge",seq_len(length(igraph::E(g))))
  #### Add nodes metadata ####
  for(x in names(vertices)){
    igraph::vertex_attr(g,x) <- as.character(
      vertices[names(igraph::V(g))][[x]]
    )
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
  #### Add edge color ####
  if(!is.null(edge_palette)){
    is_pheno_edge <- igraph::as_edgelist(g)[,1] %in% vertices[node_type=="Phenotype"]$name
    edge_color <- igraph::edge_attr(g)[grep(edge_color_var,igraph::edge_attr_names(g), value = TRUE)] |>
      data.table::as.data.table() |>
      data.table::fcoalesce()
    edge_color_dict <- stats::setNames(
      add_alpha(edge_palette(length(unique(edge_color)))),
      unique(edge_color))
    igraph::edge_attr(g,"color") <- edge_color_dict[edge_color]
    igraph::edge_attr(g,"color")[!is_pheno_edge] <- add_alpha("grey")
      igraph::edge_attr(g,"dashes") <- !is_pheno_edge
  }
  #### Add edge size ####
  if(!is.null(edge_size_var)){
    edge_size <- igraph::edge_attr(g)[grep(edge_size_var,igraph::edge_attr_names(g), value = TRUE)] |>
      data.table::as.data.table() |>
      data.table::fcoalesce()
    igraph::edge_attr(g,"width") <- edge_size
  }
  #### Format ####
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


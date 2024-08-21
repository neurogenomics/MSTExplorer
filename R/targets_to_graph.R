targets_to_graph <- function(top_targets,
                             vertex_vars,
                             group_var,
                             metadata_vars=c("hpo_id",
                                             "hpo_name",
                                             "definition",
                                             "ontLvl",
                                             "disease_name",
                                             "disease_id",
                                             "ancestor_name",
                                             "CellType",
                                             "q",
                                             "effect"
                                             # "tier_merge",
                                             # "disease_characteristic",
                                             # "gene_biotype"
                                             # grep(paste("_mean$",
                                             #            "_min$",
                                             #            "_latest$",
                                             #            "_names$",
                                             #            sep = "|"),
                                             #      names(top_targets),
                                             #      value = TRUE)
                                             ) |> unique(),
                             edge_color_var = "effect",
                             edge_size_var = edge_color_var,
                             mediator_var = "gene_symbol",
                             agg_fun=mean,
                             format="visnetwork",
                             verbose=TRUE){
  # devoptera::args2vars(targets_to_graph)
  requireNamespace("igraph")
  requireNamespace("tidygraph")
  effect <- node_type <- shape <- node <- disease_name <-
    disease_id <- hpo_name <- ancestor_name <- from <- to <- NULL;

  messager("Creating network.",v=verbose)
  #### Create vertices ####
  vertex_vars <- unique(vertex_vars[vertex_vars %in% names(top_targets)])
  shapes <- c("database","circle","box")
  if(length(vertex_vars)>3){
    shapes <- c(rep("database",length(vertex_vars)-3),shapes)
  }
  ##### Remove Phenotypes that are also ancestor #####
  ## This avoids duplicate nodes
  if("ancestor_name" %in% names(top_targets)){
    top_targets <- top_targets[hpo_name!=ancestor_name,]
  }
  if("disease_name" %in% vertex_vars){
    top_targets[,disease_name:=data.table::fcoalesce(disease_name,disease_id)]
  }
  #### Aggregate ####
  if(!is.null(agg_fun)){
    top_targets[,c(edge_color_var,edge_size_var):=list(agg_fun(get(edge_color_var)), agg_fun(get(edge_size_var))), by=vertex_vars]
  }

  #### Make vertex metadata ####
  metadata_vars <- metadata_vars[metadata_vars %in% names(top_targets)]
  lcols <- names(top_targets)[unlist(lapply(top_targets, methods::is, "list"))]
  metadata_vars <- metadata_vars[!metadata_vars %in% lcols]
  vertices <- (
    data.table::melt.data.table(
      top_targets[, unique(c(vertex_vars,group_var,metadata_vars)),with=FALSE],
      id.vars = unique(c(group_var,metadata_vars)),
      measure.vars = vertex_vars,
      variable.name = "node_type",
      value.name = "node") |>
      rev() |> unique()
  )
  vertices <- vertices[!is.na(node)]
  #### Add node shapes ####
  shape_dict <- stats::setNames(
    shapes,
    unique(vertex_vars))
  vertices[,shape:=shape_dict[node_type]]
  #### Ensure each node only appears once in the node metadata ####
  vertices <- vertices[,utils::head(.SD, 1), by = c("node")]
  #### ancestor_name is only relevant metadata for hpo_name nodes ####
  if("ancestor_name" %in% names(vertices)){
    vertices[node_type!="hpo_name",]$ancestor_name <- NA
  }
  if("definition" %in% names(vertices)){
    vertices[!node_type %in% c("hpo_name","ancestor_name")]$definition <- NA
  }
  if("ontLvl" %in% names(vertices)){
    vertices[!node_type %in% c("hpo_name","ancestor_name")]$ontLvl <- NA
  }
  vertices <- unique(vertices)
  vertices$name <- stringr::str_wrap(
    gsub("_"," ",gsub("/"," / ",vertices$node)),
    width = 10)
  #### Merge graphs ####
  if(is.character(mediator_var)){
    edges <- lapply(vertex_vars[vertex_vars!=mediator_var], function(v){
      vv <- c(v,mediator_var)
      dt <- unique(
        top_targets[,c(vv,c(edge_color_var,edge_size_var)),with=FALSE]
      )
      names(dt) <- c("from","to","color","width")
      dt
    })
  } else if(is.list(mediator_var)){
    ilist <- if(length(mediator_var)==0){
      if(length(vertex_vars)==4){
        list(c(1,2),c(2,3),c(3,4),c(2,4))
      } else if(length(vertex_vars)==5){
        list(c(1,2),c(2,3),c(3,4),c(4,5),c(3,5))
      } else if(length(vertex_vars)==6){
        list(c(1,2),c(2,3),c(3,4),c(4,5),c(5,6),c(4,6))
      }
    } else {
      mediator_var
    }
    if(length(unique(unlist(ilist)))>length(vertex_vars)){
      stp <- paste("When mediator_var is a list, mediator_var must be",
                   "equal to or less than the length of vertex_vars.")
      stop(stp)
    }
    edges <- lapply(ilist, function(il){
      vv <- vertex_vars[il]
      cols <- edge_size_var
      dt <- unique(
        top_targets[,c(vv,c(edge_color_var,edge_size_var)),
                    with=FALSE][,(cols):=lapply(.SD,mean),.SDcols=cols,by=vv]
      )
      names(dt) <- c("from","to","color","width")
      dt
    }) |> data.table::rbindlist(fill = TRUE)
  } else {
    edges <- lapply(seq(length(vertex_vars)-1), function(i){
      vv <- vertex_vars[c(i,i+1)]
      dt <- unique(
        top_targets[,c(vv,c(edge_color_var,edge_size_var)),
                    with=FALSE][,effect:=mean(effect), by=vv]
      )
      names(dt) <- c("from","to","color","width")
      dt
    }) |> data.table::rbindlist(fill = TRUE)
  }
  edges <- edges[!is.na(from) & !is.na(to)]
  #### Merge subgraphs ####
  g <- tidygraph::tbl_graph(nodes = vertices,
                            edges = edges,
                            node_key = "node")
  #### Add hoverdata ####
  g <- KGExplorer::add_hoverboxes(g = g,
                                  hoverbox_column = "title")
  #### Format ####
  if(isTRUE(format=="ggnetwork")){
    g2 <- KGExplorer::graph_to_ggnetwork(g)
    rownames(g2) <- paste0("edge",seq(nrow(g2)))
    return(g2)
  } else{
    return(g)
  }
}


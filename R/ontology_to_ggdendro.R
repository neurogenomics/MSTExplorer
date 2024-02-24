ontology_to_ggdendro <- function(ont=KGExplorer::get_ontology(
  name = "cl",
  add_ancestors = 1),
  terms=NULL #as.character(unique(results$cl_id))
  ){
  #### Get celltype dendrogram ####
  ## From CTD
  # ctd1 <- get_data("ctd_HumanCellLandscape.rds")
  # ctd2 <- get_data("ctd_DescartesHuman.rds")
  # dend_list <- EWCE:::prep_dendro(ctdIN = ctd[[1]])
  dst <- KGExplorer::ontology_to(ont = ont,
                               # terms = ids,
                               to="igraph_dist")
  terms <- intersect(terms,rownames(dst))
  # missing_terms <- terms[!terms %in% rownames(dst)]
  if(!is.null(terms)){
    dst <- stats::as.dist(dst[terms,terms])
  }

  #### Convert to hclust ####
  hc <- stats::hclust(dst)
  #### Convert to dendrogram ####
  dend <- stats::as.dendrogram(hc)
  # messager("Pruning dendrogram.")
  # dend2 <- dendextend::prune(
  #   dend,
  #   leaves=labels(dend)[!labels(dend) %in% results$cl_id])
  ddata <- ggdendro::dendro_data(dend)
  ddata$labels$id <- ddata$labels$label
  ddata$labels$label <- KGExplorer::map_ontology_terms(ont = ont,
                                                       terms = ddata$labels$label,
                                                       to = "name")
  return(ddata)
}

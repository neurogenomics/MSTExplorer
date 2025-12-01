validate_associations_mkg_dist <- function(results,
                                           kg=get_data("monarch_kg_cells.csv"),
                                           cl=get_cl(),
                                           q_threshold=0.05,
                                           alt_ids=c("CL:0000111"="CL:2000032"),
                                           metric = c("dist_nca.min_adj",
                                                      "dist_lca.min_adj",
                                                      "dist_nca.min",
                                                      "dist_lca.min")[1]
                                           ){
  from <- to <- cl_id <- dist_nca.min <- all_dist_nca.min <-
    dist_lca.min <- all_dist_lca.min <- dist_nca.min_adj <- dist_lca.min_adj <-
    hpo_id <- pct <- NULL;

  results <- map_celltype(results)
  kg <- kg[grepl("HP:",from)][from %in% unique(results$hpo_id)]
  results$cl_id <- as.character(results$cl_id)
  kg$to <- as.character(kg$to)
  #### Get distances ####
  dist_dt <- lapply(seq(nrow(kg)), function(i){
    kg_i <- kg[i,]
    res_i <- results[q<q_threshold & hpo_id==kg_i$from]
    kg_i[to %in% names(alt_ids),to:=alt_ids[to]]
    res_i[cl_id %in% names(alt_ids),cl_id:=alt_ids[cl_id]]
    # setdiff(kg_i$to,cl@terms)
    # setdiff(res_i$cl_id,cl@terms)
    messager(paste0("[i=",i,"]"),"cl_id:",kg_i$to,"& hpo_id:",kg_i$from)
    if(nrow(res_i)==0){
      return(
        data.table::data.table(
          hpo_id=kg_i$from,
          cl_id=kg_i$to
        )
      )
    } else {
      dist_nca <- simona::shortest_distances_via_NCA(
        dag = cl,
        terms = c(kg_i$to,res_i$cl_id),
        verbose=FALSE)
      all_dist_nca <- simona::shortest_distances_via_NCA(
        dag = cl,
        terms = c(kg_i$to,unique(results$cl_id)),
        verbose=FALSE)
      dist_lca <- simona::longest_distances_via_LCA(
        dag = cl,
        terms = c(kg_i$to,res_i$cl_id),
        verbose=FALSE)
      all_dist_lca <- simona::longest_distances_via_LCA(
        dag = cl,
        terms = c(kg_i$to,unique(results$cl_id)),
        verbose=FALSE)
      data.table::data.table(
        hpo_id=kg_i$from,
        cl_id=kg_i$to,
        dist_nca.min=min(dist_nca[kg_i$to, res_i$cl_id]),
        dist_nca.max=max(dist_nca[kg_i$to, res_i$cl_id]),
        dist_nca.median=stats::median(dist_nca[kg_i$to, res_i$cl_id]),

        dist_lca.min=min(dist_lca[kg_i$to, res_i$cl_id]),
        dist_lca.max=max(dist_lca[kg_i$to, res_i$cl_id]),
        dist_lca.median=stats::median(dist_lca[kg_i$to, res_i$cl_id]),

        all_dist_nca.min=min(all_dist_nca[kg_i$to, unique(results$cl_id)], na.rm = TRUE),
        all_dist_lca.min=min(all_dist_lca[kg_i$to, unique(results$cl_id)], na.rm = TRUE)
      )[,dist_nca.min_adj:=((dist_nca.min+1)/(all_dist_nca.min+1))-1][,dist_lca.min_adj:=((dist_lca.min+1)/(all_dist_lca.min+1))-1]
      ## ^ Adjust for the fact that the CellTypeDataset references don't
      ## contain every possible cell type. This makes it a more fair comparison.
    }
  })|> data.table::rbindlist(fill=TRUE, idcol = "i")
  #### Aggregate data by degree dist ####
  plot_dat <- lapply(seq(0,max(dist_dt[[metric]], na.rm = TRUE)),
                     function(d){
    data.table::data.table(
      dist=d,
      pct=nrow(dist_dt[get(metric)<=d,])/nrow(dist_dt)*100
    )
  }) |> data.table::rbindlist()
  #### Create plot ####
  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x=as.factor(dist),
                                         y=pct,
                                         label=paste0(round(pct,1),"%"),
                                         group=1)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_label(vjust=-.25) +
    ggplot2::labs(
      x="Ontological distance between cell types",
      y="Recall of known associations (%)"
    ) +
    ggplot2::lims(y=c(0,100)) +
    ggplot2::theme_minimal()
  #### Plot ####
  return(
    list(
      data=plot_dat,
      plot=p
    )
  )
}

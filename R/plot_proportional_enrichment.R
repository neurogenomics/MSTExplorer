plot_proportional_enrichment <- function (results=load_example_results(),
                                           target_branches = get_target_branches(),
                                           target_celltypes = get_target_celltypes(target_branches=target_branches),
                                           celltype_col="cl_id",
                                           hpo=HPOExplorer::get_hpo(),
                                           color_map = KGExplorer::map_colors(
                                             data.table::data.table(
                                               ancestor=names(target_branches)
                                               ),
                                             preferred_palettes = "tol"
                                             )[[1]],
                                           legend.position = "none",
                                           y_limits=NULL,
                                           scales="free_y",
                                           y_var="pct",
                                           y_lab=if(y_var=="enrichment") {
                                             "On-target enrichment"
                                           } else {
                                             "On-target associations (%)"
                                           },
                                           nrow=1,
                                           show_plot=TRUE) {
 #  res <- rols::OlsSearch(q = "nervous system",
 #                         ontology = "cl")
 #  res <- rols::olsSearch(res)
 # View( as(res, "data.frame"))
 #  as(res, "Terms")
 #
 #  rols::descendants(id = "UBERON:0016879")
 #
  # cl_dt <- KGExplorer::ontology_to(ont = cl,
  #                                  to = "data.table",
  #                                  what="nodes",
  #                                  prefox='')
  # results2 <- merge(results,
  #       cl_dt[,c("subject","subject_ancestor_name")][,.SD[1], by="subject"],
  #       all.x = TRUE,
  #       by.x="cl_id",
  #       by.y="subject")
  requireNamespace("ggplot2")
  ancestor_name <- pct <- pct_min <- pct_max <- enrichment <- pct_celltype <-
    enrichment_mean <- NULL;

  results <- HPOExplorer::add_hpo_name(results)
  results <- HPOExplorer::add_ancestor(results,
                                       hpo = hpo)
  results[,minus_log_p:=round(-log10(p +1e-06), digits = 0)]
  #### Add Cell Ontology labels ####
  # target_cells <- lapply(target_cells,EWCE::fix_celltype_names)
  results <- map_celltype(results)
  ### Repeat checks over each HPO (branch) ####
  n_all_celltypes <- length(unique(results[[celltype_col]]))
  proportional_results <- lapply(stats::setNames(names(target_celltypes),
                                                 names(target_celltypes)),
                                 function(target_ancestor){
    tct <- target_celltypes[[target_ancestor]]
    tct <- intersect(tct, unique(results[[celltype_col]]))
    d <- results[ancestor_name == target_ancestor]
    results[,list(
      pct=sum(ancestor_name==target_ancestor & get(celltype_col) %in% tct)/
          sum(get(celltype_col) %in% tct)*100,
      n_celltype=length(tct),
      pct_celltype=length(tct)/n_all_celltypes*100,
      total_per_celltype=sum(get(celltype_col) %in% tct)
      # celltype_intersect=paste(tct,collapse = ";")
      ),
      by=c("minus_log_p")][,celltype_label:=paste0(
        paste(target_branches[[target_ancestor]],collapse = " /\n"),
        " (n=",length(unique(tct)),")")]
  }) |> data.table::rbindlist(idcol = "branch")

  #### Add newline to long branch names #####
  proportional_results[,branch:=gsub("Abnormality of the","Abnormality of the\n",branch)]
  names(color_map) <- gsub("Abnormality of the","Abnormality of the\n",names(color_map))
  #### Compute enrichment (observed/expected) ####
  proportional_results[,enrichment:=(pct/pct_celltype)]
  #### Compute summary stats #####
  ss <- proportional_results[,list(pct_max=max(pct),
                                   pct_min=min(pct),
                                   enrichment_mean=mean(enrichment)),
                       by="branch"][,list(pct_min=max(pct_min),
                                          pct_max=max(pct_max),
                                          pct_max_mean=mean(pct_max),
                                          pct_max_sd=stats::sd(pct_max),
                                          enrichment_mean=mean(enrichment_mean)
                                          )]
  messager("Proportional enrichment summary stats:")
  messager(" - pct_min:",round(ss$pct_min,2))
  messager(" - pct_max:",round(ss$pct_max,2))
  messager(" - pct_max_mean:",round(ss$pct_max_mean,2))
  messager(" - pct_max_sd:",round(ss$pct_max_sd,2))
  messager(" - enrichment_mean:",round(ss$enrichment_mean,2))
  #### Plot ####
  baseline_dat <- proportional_results[,.SD[1],by=c("branch")]
  baseline_dat[,baseline:=ifelse(y_var=="enrichment",1,pct_celltype)]
  prop_plot <- ggplot2::ggplot(proportional_results,
                               ggplot2::aes(x = minus_log_p,
                                            y = !!ggplot2::sym(y_var),
                          color = branch),
                      na.rm = TRUE) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(nrow = nrow,
                        scales = scales,
                        facets=c("branch","celltype_label")) +
    ggplot2::labs(x=bquote(-log[10](p)),
                  y=y_lab) +
    ggplot2::scale_y_continuous(limits = y_limits) +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::geom_hline(data=baseline_dat,
                        mapping = ggplot2::aes(yintercept = baseline),
                        linetype="dashed",
                        alpha=.5) +
    ggplot2::geom_text(data = baseline_dat,
              mapping = ggplot2::aes(y = baseline,
                                     x=6,
                                     label=paste0(
                                       "baseline ",
                                                  "(",n_celltype,"/",n_all_celltypes,")"),
                            ),
              hjust=1,
              vjust = -0.75,
              size=3,
              color="black",
              alpha=.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill="transparent"),
          legend.position = legend.position)
  if(isTRUE(show_plot)) methods::show(prop_plot)
  #### Return ####
  return(list(plot=prop_plot,
              data=proportional_results))
}

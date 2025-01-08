#' Prioritise targets: multi-heatmap
#'
#' @param top_targets A results table after it has been annotated with
#' \link[MSTExplorer]{add_disease} and \link[MSTExplorer]{add_symptom_results}.
#' @param prioritise_targets_network_out Output of
#' \link[MSTExplorer]{prioritise_targets_network}.
#' @param gene_order The order in which to show genes in heatmap x-axis.
#' @param gencc_extra Extra rows to add to the diseases GenCC data.
#' @param size Plot text size.
#'
#' @inheritParams plot_
#' @inheritParams prioritise_targets
#' @inheritParams patchwork::wrap_plots
#' @inheritParams stringr::str_wrap
#' @export
#' @examples
#' top_targets <- MSTExplorer::example_targets$top_targets[1:10]
#' prioritise_targets_network_out <- prioritise_targets_network(
#'   top_targets = top_targets)
#'
#' out <- prioritise_targets_multiheatmap(
#'   top_targets = top_targets,
#'   prioritise_targets_network_out = prioritise_targets_network_out,
#'   ctd_list=load_example_ctd("ctd_DescartesHuman.rds", multi_dataset = TRUE))
prioritise_targets_multiheatmap <- function(top_targets,
                                            prioritise_targets_network_out,
                                            ctd_list=load_example_ctd(
                                              c("ctd_DescartesHuman.rds",
                                                "ctd_HumanCellLandscape.rds"),
                                              multi_dataset = TRUE),
                                            gene_order=NULL,
                                            hpo=HPOExplorer::get_hpo(),
                                            gencc_extra=list(),
                                            # gencc_extra=list("disease_id"="OMIM:614372",
                                            #                  "disease_name"="Mannose-binding lectin (MBL) deficiency",
                                            #                  "gene_symbol"="MBL2",
                                            #                  "evidence_score_sum"=5),
                                            size=6,
                                            width=50,
                                            heights = NULL,
                                            # heights = c(7/7,1/7,4/7)
                                            show_plot=TRUE
                                            ){
  node_type <- gene_symbol <- disease_id <- specificity <- evidence_score_sum <-
    disease_name <- hpo_name <- mean_specificity <- cl_name <- NULL;

  ctd_list <- ctd_list[unique(top_targets$ctd)]
  if(length(ctd_list)==0){
    stopper("No matching CTDs found in both `ctd_list` and `top_targets$ctd`.")
  }
  top_targets <- map_celltype(top_targets)
  graph_dat <- KGExplorer::graph_to_dt(prioritise_targets_network_out$data,
                                       what = "nodes")
  phenotype <- unique(graph_dat$hpo_name)
  phenotype_id  <- unique(graph_dat$hpo_id)
  if(is.null(gene_order)){
    gene_order <- graph_dat[node_type=="gene_symbol"]$node
  }
  gencc_disease <- (
    KGExplorer::get_gencc()
  )[disease_id %in% unique(graph_dat$disease_id)] |>
    HPOExplorer::map_disease()|>
    ## Add data that's now missing from GenCC
    rbind(
      gencc_extra,
      fill=TRUE
    )|>
    dplyr::mutate(gene_symbol=factor(gene_symbol, gene_order, ordered=TRUE))

  gencc_phenotype <- (
    HPOExplorer::hpo_to_matrix(terms = phenotype_id,
                               as_matrix = FALSE) |>
      data.table::melt.data.table(id.vars = "gene_symbol",
                                  variable.name = "hpo_id",
                                  value.name = "evidence_score")
  )|>
    HPOExplorer::add_hpo_name(hpo = hpo)|>
    dplyr::filter(gene_symbol %in% gene_order)|>
    dplyr::mutate(gene_symbol=factor(gene_symbol, gene_order, ordered=TRUE))



  get_ctd_matrix <- function(ctd,
                             lvl,
                             ctd_name,
                             metric="specificity",
                             rows=NULL,
                             cols=NULL){
    X <- ctd[[lvl]][[metric]]
    if(is.null(rows)) {
      rows <- rownames(X)
    } else{
      rows <- intersect(rows, rownames(X))
    }
    if(is.null(cols)) {
      cols <- colnames(X)
    } else {
      cols <- intersect(cols, colnames(X))
    }
    (
      X[rows,cols] |>
        as.matrix() |>
        data.table::as.data.table(keep.rownames="gene_symbol") |>
        data.table::melt.data.table(id.vars = "gene_symbol",
                                    variable.name = "CellType",
                                    value.name = metric))[,ctd:=eval(ctd_name)]|>
        map_celltype()
  }


  top_targets_sub <- top_targets[hpo_id %in% graph_dat$hpo_id &
                                   disease_id %in% graph_dat$disease_id &
                                   gene_symbol %in% gene_order]
  ctd_rni <- (
    lapply(names(ctd_list), function(x){
      lvl <- map_ctd_levels(top_targets_sub)[x]$annotLevel
      get_ctd_matrix(ctd=ctd_list[[x]],
                     ctd_name=x,
                     rows=gene_order,
                     cols=top_targets_sub[ctd==eval(x)]$CellType,
                     lvl=lvl)
    }) |> data.table::rbindlist()
  )[,list(mean_specificity=mean(specificity)), by=c("cl_name","gene_symbol")]|>
    dplyr::mutate(gene_symbol=factor(gene_symbol, gene_order, ordered=TRUE))


  # top_targets_sub <- MSTExplorer:::add_driver_genes(
  #   results = top_targets[hpo_id %in% graph_dat$hpo_id &
  #                           disease_id %in% graph_dat$disease_id &
  #                           gene_symbol %in% graph_dat[node_type=="gene_symbol"]$name],
  #   ctd_list = ctd_list,
  #   metric = "specificity")
  # top_targets_sub[,mean_specificity:=mean(specificity), by=c("cl_name","gene_symbol")]|>
  #   dplyr::mutate(gene_symbol=factor(gene_symbol, gene_order, ordered=TRUE))


  heatmap_disease <- ggplot2::ggplot(gencc_disease,
                                     ggplot2::aes(x=gene_symbol,
                                                  y=stringr::str_wrap(disease_name,width),
                                                  fill=evidence_score_sum)) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic() +
    ggplot2::labs(x=NULL,y="Disease",fill="Evidence\nscore") +
    # ggplot2::scale_fill_viridis_c() +
    ggplot2::scale_fill_distiller(palette = "Blues",
                                  limits=c(0,max(gencc_disease$evidence_score_sum, na.rm = TRUE)),
                                  n.breaks=4,
                                  na.value = "transparent",
                                  direction = 1) +
    # ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = barheight)) +
    ggplot2::theme(text = ggplot2::element_text(size=size),
                   legend.title = ggplot2::element_text(size=size),
                   legend.text = ggplot2::element_text(size=size),
                   axis.text.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin())

  heatmap_phenotype <- ggplot2::ggplot(gencc_phenotype,
                                       ggplot2::aes(x=gene_symbol,
                                                    y=stringr::str_wrap(hpo_name,width),
                                                    fill=evidence_score)
  ) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic() +
    ggplot2::labs(x=NULL,y="Phenotype",fill="Evidence\nscore") +
    # ggplot2::scale_fill_viridis_c() +
    ggplot2::scale_fill_distiller(palette = "Purples",
                                  limits=c(0,max(gencc_phenotype$evidence_score)),
                                  na.value = "transparent",
                                  n.breaks=4,
                                  direction = 1) +
    # ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = barheight)) +
    ggplot2::theme(text = ggplot2::element_text(size=size),
                   legend.title = ggplot2::element_text(size=size),
                   legend.text = ggplot2::element_text(size=size),
                   axis.text.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin())

  heatmap_celltype <- ggplot2::ggplot(ctd_rni,
                                      ggplot2::aes(x=gene_symbol,
                                                   y=stringr::str_wrap(cl_name,width),
                                                   fill=mean_specificity)) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic() +
    ggplot2::labs(x="Gene",y="Cell type",fill="Expression\nspecificity") +
    # ggplot2::scale_fill_viridis_c() +
    ggplot2::scale_fill_distiller(palette = "YlOrBr",
                                  values=c(.Machine$double.xmin,1),
                                  limits=c(0,1),
                                  n.breaks=4,
                                  na.value = "transparent",
                                  direction = 1) +
    # ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = barheight)) +
    ggplot2::theme(text = ggplot2::element_text(size=size),
                   legend.title = ggplot2::element_text(size=size),
                   legend.text = ggplot2::element_text(size=size),
                   plot.margin = ggplot2::margin())
  ## Merge plots
  if(is.null(heights)){
    heights <- c(
      length(unique(gencc_disease$disease_name)),
      length(unique(gencc_phenotype$hpo_name)),
      length(unique(ctd_rni$cl_name))
    )
    heights <- heights / sum(heights)
  }
  ggp <- patchwork::wrap_plots(heatmap_disease,
                               heatmap_phenotype,
                               heatmap_celltype,
                               heights = heights,
                               ncol = 1)
  if(show_plot) methods::show(ggp)
  return(list(
    data=list(
      "Disease"=gencc_disease,
      "Phenotype"=gencc_phenotype,
      "Cell type"=ctd_rni
    ),
    plot=ggp
  ))
}

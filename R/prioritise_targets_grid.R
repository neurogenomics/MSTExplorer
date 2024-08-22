#' Prioritise targets grid
#'
#' Plot the output of \link[MSTExplorer]{prioritise_targets}
#' as a grid.
#' @param top_targets output of \link[MSTExplorer]{prioritise_targets}.
#' @param species_list Species to include in orthologous genes grid.
#' @param keep_severity_class Phenotypes to keep based on severity classes.
#' @param keep_physical_malformations Phenotypes to keep based on physical
#' malformation frequency (0=never, 1=rarely, 2=often, 3=always).
#' @param widths1 Proportional widths of severity annotation heatmap
#'  (left subplot) and the top targets/orthologous genes grid (right subplot).
#' @param widths2 Proportional widths of the top targets/orthologous genes grid.
#' @inheritParams plot_
#' @inheritParams prioritise_targets
#' @inheritParams HPOExplorer::plot_top_phenos
#' @inheritDotParams HPOExplorer::plot_top_phenos
#' @inheritParams patchwork::plot_layout
#' @inheritParams stringr::str_trunc
#' @returns A named list with data and patchwork object.
#'
#' @export
#' @importFrom orthogene convert_orthologs
#' @examples
#' top_targets <- MSTExplorer::example_targets$top_targets
#' out <- prioritise_targets_grid(top_targets = top_targets)
prioritise_targets_grid <- function(top_targets,
                                    res_class = HPOExplorer::gpt_annot_class(),
                                    n_per_class = 10,
                                    keep_severity_class=c("profound","severe"),
                                    keep_physical_malformations=NULL,
                                    species_list = c("Homo sapiens",
                                                     "Macaca mulatta",
                                                     "Mus musculus",
                                                     "Danio rerio",
                                                     "Drosophila melanogaster",
                                                     "Caenorhabditis elegans"),
                                    legend.position='left',
                                    keep_ont_levels=NULL,
                                    width=70,
                                    widths1=c(1,4),
                                    widths2=c(1,.8),
                                    show_plot=TRUE,
                                    ...
                                    ){
  hpo_name <- severity_class <- physical_malformations <- value <- variable <-
    p <- species <- NULL;

  top_targets <- top_targets|> data.table::copy()
  top_targets <- map_celltype(top_targets)
  top_targets <- HPOExplorer::add_disease(top_targets,
                                          add_descriptions = TRUE)
  #### Filter severity classes ####
  if(!is.null(keep_severity_class)){
    res_class <- res_class[severity_class %in% keep_severity_class]
  }
  if(!is.null(keep_physical_malformations)){
    res_class <- res_class[physical_malformations %in% keep_physical_malformations]
  }
  #### Filter by intersect between severity metdata and top_targets ####
  select_phenos <- intersect(res_class$hpo_name,
                             unique(top_targets$hpo_name))
  res_class <- res_class[hpo_name %in% select_phenos]
  if(nrow(res_class)==0){
    stopper("No phenotypes remaining after filtering.")
  }
  # prioritise_targets_out$top_targets
  plot_top_phenos_out <- HPOExplorer::plot_top_phenos(
    res_class=res_class,
    keep_ont_levels=keep_ont_levels,
    legend.position=legend.position,
    nrow=2,
    n_per_class = n_per_class,
    ...)

  tile_plot <- function(dt1,
                        dt2,
                        subtitle=c("Top targets","Orthologous genes"),
                        label_cols=c("Disease"="disease_name",
                                     "Cell type"="cl_name"),
                        x_pos=c(1,2),
                        labels=names(label_cols),
                        xtext=TRUE,
                        width=70,
                        widths=c(1,.8),
                        size=3){
    if(!is.null(labels)) labels <- stats::setNames(labels,as.character(x_pos))
    dat <- merge(dt1[,-c("hpo_id")],
                 dt2, by="hpo_name", sort=FALSE)|>
      data.table::melt.data.table(
        id.vars=c("hpo_id","hpo_name","severity_class","p"),
        measure.vars=unique(c(label_cols,"gene_symbol")))
    dat[,hpo_name:=factor(hpo_name,
                          levels=rev(levels(dt1$hpo_name)), ordered=TRUE)]
    dat[startsWith(value,"GRANULOMATOUS"), value:=stringr::str_to_title(value)]
    dat <- dat[,.SD[1], by=c("hpo_name","variable")]
    dat[,x:=stats::setNames(x_pos, label_cols)[variable]]
    dat[p==0,p:=.Machine$double.xmin]

    ggtile <- ggplot2::ggplot(dat[variable %in% label_cols],
                              ggplot2::aes(
                                x=x, y=hpo_name,
                                fill=-log1p(p),
                                label=stringr::str_trunc(value,width=width),
                                hpo_id=hpo_id)) +
      ggplot2::scale_y_discrete(drop=TRUE) +
      ggplot2::geom_text(hjust=0, size=size) +
      ggplot2::theme_classic() +
      ggplot2::scale_x_continuous(labels=labels, breaks=x_pos,
                                  limits=c(min(x_pos), max(x_pos)*1.2)) +
      ggplot2::facet_grid(severity_class~., scales = "free_y") +
      ggplot2::labs(x=NULL, y=NULL, subtitle=subtitle[1]) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     strip.text=ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_line(),
                     panel.grid.minor.y = ggplot2::element_line() )
    if(isFALSE(xtext)) {
      ggtile <- ggtile + ggplot2::theme(axis.text.x = ggplot2::element_blank())
    }
    #### Gene ortholog plot ####
    orthology <- lapply(stats::setNames(species_list,species_list), function(s){
      if(s=="Homo sapiens") return( dat[variable=="gene_symbol",
                                        list(input_gene=value,
                                             ortholog_gene=value)] )
      orthogene::convert_orthologs(dat[variable=="gene_symbol"],
                                   gene_input="value",
                                   input_species="human",
                                   output_species=s,
                                   # non121_strategy="kbs",
                                   gene_output="column",
                                   method="homologene")
    })|> data.table::rbindlist(idcol="species")
    orthology[,species:=factor(species,levels=species_list, ordered=TRUE)]
    dt <- merge(dat[variable=="gene_symbol"],
                orthology,
                by.x="value",
                by.y="input_gene",
                allow.cartesian=TRUE,
                all.x=TRUE)
    gg <- ggplot2::ggplot(dt, ggplot2::aes(x=species, y=hpo_name,
                                           label=ortholog_gene)) +
      ggplot2::scale_x_discrete(labels=stringr::str_wrap(species_list,10),
                                drop=FALSE) +
      ggplot2::geom_tile(color="grey20", fill="white") +
      ggplot2::geom_text(color="black", size=size, check_overlap = TRUE) +
      ggplot2::theme_minimal() +
      ggplot2::facet_grid(severity_class~., scales = "free_y") +
      ggplot2::theme(#axis.text.x=ggplot2::element_text(angle=45, hjust=1),
        axis.text.y=ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        strip.text=ggplot2::element_blank()) +
      ggplot2::labs(x=NULL,y=NULL, subtitle=subtitle[2]) +
      ggplot2::scale_fill_viridis_c()
    if(isFALSE(xtext)) {
      gg <- gg + ggplot2::theme(axis.text.x=ggplot2::element_blank())
    }
    ggtile <- (ggtile | gg) + patchwork::plot_layout(widths=widths)

    return(list(
      data=list(top_targets=dat,
                orthologous_genes=dt),
      plot=ggtile
    ))
  }

  tp1 <- tile_plot(dt1=plot_top_phenos_out$data$congenital,
                   dt2=top_targets,
                   width=width,
                   widths=widths2,
                   xtext=FALSE)
  tp2 <- tile_plot(dt1=plot_top_phenos_out$data$noncongenital,
                   dt2=top_targets,
                   width=width,
                   widths=widths2,
                   subtitle=c("",""))
  gg_labels <- tp1$plot/tp2$plot
  gg_top_phenos <- (plot_top_phenos_out$plot | gg_labels) +
    patchwork::plot_layout(widths=widths1)
  #### Show plot ####
  if(show_plot) methods::show(gg_top_phenos)
  #### Return ####
  return(
    list(data=list(congenital=tp1$data,
                   noncongenital=tp2$data),
         plot=gg_top_phenos)
  )
}

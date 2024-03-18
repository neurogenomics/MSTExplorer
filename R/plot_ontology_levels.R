#' Plot ontology levels
#'
#' Generate plots comparing the ontology level of each HPO phenotype and
#' several other metrics.
#' @param group_vars Compute the mean value for each \code{x_vars}
#'  grouped by the \code{group_vars}.
#'  This can be less accurate for some metrics but helps to drastically
#'  reduce computational load.
#' @param p2g Phenotype to gene data.
#' @param x_vars Variables to plot on the x-axis of each subplot.
#' @param min_value Minimum value for the \code{specificity} metric.
#' @param n.breaks Passed to \link[ggplot2]{scale_fill_viridis_c}.
#' @param log_vars Logical vector indicating which variables to log-transform.
#' @param sig_vars Logical vector indicating which variables to only plot
#' for significant results.
#' @inheritParams plot_
#' @inheritParams ggpubr::stat_cor
#' @inheritParams prioritise_targets
#' @inheritParams ggplot2::geom_boxplot
#' @inheritParams ggstatsplot::ggscatterstats
#' @returns A named list containing the data and the plot.
#'
#' @export
#' @import HPOExplorer
#' @examples
#' out <- plot_ontology_levels()
plot_ontology_levels <- function(results = load_example_results(),
                         p2g = HPOExplorer::load_phenotype_to_genes(),
                         ctd_list = load_example_ctd(
                           file = paste0("ctd_",unique(results$ctd),".rds"),
                           multi_dataset = TRUE
                         ),
                         x_vars = c("genes",
                                    "cell types",
                                    "estimate",
                                    "mean_specificity"),
                         log_vars=x_vars %in% c("estimate","statistic",
                                                "F","ges",
                                                "effect"),
                         sig_vars=x_vars %in% c("estimate","statistic",
                                                "F","ges",
                                                "effect",
                                                "cell types",
                                                "mean_specificity"),
                         group_vars=c("hpo_id",
                                      "ontLvl",
                                      "ctd"),
                         q_threshold=0.05,
                         min_value=NULL,
                         label.x.npc = .05,
                         label.y.npc = .5,
                         n.breaks = 4,
                         notch = FALSE,
                         nrow = 2,
                         show_plot=TRUE,
                         save_path=NULL,
                         height=7,
                         width=length(x_vars)*5.75,
                         smooth.line.args=list(method = "loess",
                                               se = FALSE)
                         ){

  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  requireNamespace("ggpubr")
  requireNamespace("ggstatsplot")
  requireNamespace("gginnards")
  gene_symbol <- CellType <- specificity <- `cell types` <- NULL;

  #### Check input variables
  ## log_vars
  if(length(log_vars)!=length(x_vars)){
    stopper("log_vars must be the same length as x_vars")
  }
  ## sig_vars
  if(length(sig_vars)!=length(x_vars)){
    stopper("sig_vars must be the same length as x_vars")
  }
  #### Preprocess data ####
  results <- HPOExplorer::add_ont_lvl(results)
  q_tmp <- if(sig_vars[which(x_vars=="cell types")]) q_threshold else 1
  results[,`cell types`:=length(unique(CellType[q<q_tmp])),
          by=group_vars]
  #### Add specificity of DRIVER GENES ####
  if(sum(grepl("specificity$",x_vars))>0){
    driver_genes <- add_driver_genes(results[q<q_threshold],
                                     ctd_list = ctd_list,
                                     metric="specificity",
                                     min_value=min_value,
                                     allow.cartesian=TRUE)
    driver_dt <- driver_genes[,list(
      driver_genes=length(unique(gene_symbol)),
      max_specificity=max(specificity),
      mean_specificity=mean(specificity),
      min_specificity=min(specificity)),
      by=c("hpo_id","CellType")]
    results <- merge(results,
                     driver_dt,
                     all.x=TRUE,
                     by=c("hpo_id","CellType"))
  }
  #### Add total gene count ####
  pcount <- p2g[,list(genes=length(unique(gene_symbol))), by="hpo_id"]
  r2 <- merge(results,
              pcount,
              by="hpo_id")
  group_vars <- group_vars[group_vars %in% names(r2)]
  #### iterator plot function
  plot_func <- function(x_var="genes",
                        y_var="ontLvl",
                        geom="boxplot",
                        method = "loess",
                        span = 0.75,
                        direction = 1,
                        trans=NULL,
                        reduce_fun=mean,
                        type ="nonparametric",
                        ...){
    # devoptera::args2vars(plot_func)
    y_lab <- x_var
    title <- paste("Phenotype level vs.",x_var)
    #### Filter to only sig vars ####
    if(sig_vars[which(x_vars==x_var)]){
      r2 <- data.table::copy(r2[q<q_threshold,])
      group_vars <- union(group_vars,"CellType")
      title <- paste(title,"(FDR<0.05)")
    }
    #### Log vars ####
    if(log_vars[which(x_vars==x_var)]){
      trans <- "log10"
      ## ggstatsplot doesn't understand how to handle expressions...
      y_lab <- paste0("log10(",y_lab,")")
      # y_lab <- substitute(expression(log[10](y_lab)),
      #                     list(y_lab=y_lab))
    }
    if(!is.null(reduce_fun)){
      dat <- r2[!is.na(get(x_var)),
                list(reduce_fun(get(x_var)))|>
                  `names<-`(x_var), by=group_vars]
    } else {
      dat <- r2[!is.na(get(x_var))]
    }
    #### Compute mean value per ont level for color #####
    dat[,mean:=mean(get(x_var), na.rm=TRUE),by=y_var]
    ####
    if(geom=="ggscatterstats"){
      p <- ggstatsplot::ggscatterstats(data = dat,
                                       x = !!ggplot2::sym(y_var),
                                       y = !!ggplot2::sym(x_var),
                                       marginal=FALSE,
                                       bf.message=FALSE,
                                       # k = 1L,
                                       type =type,
                                       point.args = list(alpha=.1),
                                       title=title,
                                       smooth.line.args=smooth.line.args) +
        ggplot2::geom_boxplot(orientation = "x",
                              ggplot2::aes(group= !!ggplot2::sym(y_var),
                                           fill=mean),
                              position = "dodge",
                              outlier.alpha = 0) +
        ggplot2::scale_fill_viridis_c(option = "plasma",
                                      n.breaks = n.breaks,
                                      direction = direction) +
        ggplot2::coord_flip() +
        ggplot2::scale_x_reverse() +
        ggplot2::labs(x="Phenotype ontology level",
                      y=y_lab,
                      fill="Per-level mean"
                      ) +
        ggplot2::theme(legend.position = "bottom",
                       legend.key.width=ggplot2::unit(1,"cm"),
                       plot.subtitle = ggplot2::element_text(size=10)
      )
      p <- gginnards::move_layers(p,idx = 2, position = "top")
      return(p)
    }
    #### Switch plots ####
    gp <- ggplot2::ggplot(dat,
                          na.rm = TRUE,
                          ggplot2::aes(x=!!ggplot2::sym(y_var),
                                       y=!!ggplot2::sym(x_var),
                                       fill=mean
                                       # color=ctd
                                       )
                          ) +
      ggplot2::scale_fill_viridis_c(option = "plasma",
                                    n.breaks = n.breaks,
                                    direction = direction)
    # if(x_var=="effect") trans <- "log"
    if(!is.null(trans)){
      gp <- gp + ggplot2::scale_x_continuous(trans = trans)
    }
    if(geom=="hex"){
      gp <- gp + ggplot2::geom_hex(...)
    } else if(geom=="violin"){
      gp <- gp + ggplot2::geom_violin(orientation = "x",
                                      ggplot2::aes(fill=ontLvl),
                             ...)
    } else if(geom=="boxplot"){
      gp <- gp + ggplot2::geom_boxplot(orientation = "x",
                                     ggplot2::aes(group=ontLvl),
                                     position = "dodge",
                                     outlier.alpha = 0.25,
                              ...)
    } else if (geom=="jitter") {
      gp <- gp + ggplot2::geom_jitter(width = 0,
                                      alpha=.25,
                            ...)
    }
    gp <- gp +
      ggplot2::geom_smooth(method = method,
                           span = span,
                           se = FALSE) +
      ggpubr::stat_cor(method = "pearson",
                       label.x.npc = label.x.npc,
                       label.y.npc = label.y.npc) +
      ggplot2::coord_flip() +
      ggplot2::scale_x_reverse() +
      ggplot2::labs(x="Phenotype ontology level",
                    y=y_lab,
                    title=title) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top",
                     legend.key.width=ggplot2::unit(1,"cm"))

    return(gp)
  }
  #### Create plots ####
  plts2 <- lapply(stats::setNames(x_vars,x_vars),
                  plot_func,
                  y_var="ontLvl",
                  geom="ggscatterstats",
                  type="parametric",
                  method = "loess",
                  span = .5,
                  direction = 1,
                  notch=notch)
  #### Merge plots ####
  plts2 <- patchwork::wrap_plots(plts2, nrow = nrow) +
    patchwork::plot_annotation(tag_levels = letters[seq_len(length(plts2))]) +
    patchwork::plot_layout(axis_titles = "collect",
                           axes="collect")
  #### Extract ggstatsplot results ####
  data_stats <- get_ggstatsplot_stats(plts2)
  #### Show plot ####
  if(show_plot) methods::show(plts2)
  #### Save plot ####
  if(!is.null(save_path)){
    KGExplorer::plot_save(plt = plts2,
                          save_path = save_path,
                          height = height,
                          width = width)
  }
  return(list(data=r2,
              data_stats=data_stats,
              plot=plts2))
}


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
#' @inheritParams ggpubr::stat_cor
#' @returns A named list containing the data and the plot.
#' @inheritParams prioritise_targets
#'
#' @export
#' @import HPOExplorer
#' @examples
#' plts <- plot_ont_lvl(label.x.npc = .05)
plot_ont_lvl <- function(results = load_example_results(multi_dataset = TRUE),
                         p2g = HPOExplorer::load_phenotype_to_genes(),
                         x_vars = c("genes",
                                    "celltypes",
                                    "log_fold_change"),
                         group_vars=c("hpo_id",
                                      "ontLvl",
                                      "ctd"),
                         label.x.npc = .05,
                         label.y.npc = .5,
                         notch = FALSE
                         ){

  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  requireNamespace("ggpubr")
  gene_symbol <- celltypes <- CellType <-
    log_fold_change <- fold_change <- NULL;

  results <- HPOExplorer::add_ont_lvl(results)
  results[,celltypes:=length(unique(CellType[q<0.05])), by=group_vars]
  pcount <- p2g[,list(genes=length(unique(gene_symbol))), by="hpo_id"]

  r2 <- merge(results,
              pcount,
              by="hpo_id")
  r2[,log_fold_change:=log(abs(fold_change))]
  group_vars <- group_vars[group_vars %in% names(r2)]
  # r2=r2[1:50000]


  plot_func <- function(x_var="genes",
                        y_var="ontLvl",
                        geom="boxplot",
                        method = "loess",
                        direction = 1,
                        trans=NULL,
                        reduce_fun=mean,
                        ...){
    #### Only include logFC values with FDR<5% ####
    if(x_var=="log_fold_change"){
      r2 <- data.table::copy(r2[q<0.05,])
    }
    if(!is.null(reduce_fun)){
      dat <- r2[,list(reduce_fun(get(x_var)))|> `names<-`(x_var), by=group_vars]
    } else {
      dat <- r2
    }
    #### Compute mean value per ont level for color #####
    dat[,mean:=mean(get(x_var)),by=y_var]
    gp <- ggplot2::ggplot(dat,
                          na.rm = TRUE,
                          ggplot2::aes(x=!!ggplot2::sym(y_var),
                                       y=!!ggplot2::sym(x_var),
                                       fill=mean
                                       # color=ctd
                                       )
                          ) +
      ggplot2::scale_fill_viridis_c(option = "plasma",
                                    direction = direction)
    # if(x_var=="fold_change") trans <- "log"
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
    } else {
      gp <- gp + ggplot2::geom_jitter(width = 0,
                                      alpha=.25,
                            ...)
    }
    gp <- gp +
      ggplot2::geom_smooth(method = method,
                           se = FALSE) +
      ggpubr::stat_cor(method = "pearson",
                       label.x.npc = label.x.npc,
                       label.y.npc = label.y.npc) +
      ggplot2::theme_bw() +
      ggplot2::coord_flip() +
      ggplot2::scale_x_reverse()

    return(gp)
  }
  # plts1 <- lapply(c("log(genes)","log(celltypes)"), plt) |>
  #   patchwork::wrap_plots(ncol = 1)
  plts2 <- lapply(stats::setNames(x_vars,x_vars),
                  plot_func,
                  y_var="ontLvl",
                  geom="boxplot",
                  direction = 1,
                  notch=notch)
  plts2 <- patchwork::wrap_plots(plts2, ncol = 1)
  return(list(data=r2,
              plot=plts2))
}


#' Plot ontology levels
#'
#' Generate plots comparing the ontology level of each HPO phenotype and
#' several other metrics.
#' @param p2g Phenotype to gene data.
#' @param x_vars Variables to plot on the x-axis of each subplot.
#' @returns A named list containing the data and the plot.
#' @inheritParams prioritise_targets
#'
#' @export
#' @import HPOExplorer
#' @examples
#' plts <- plot_ont_lvl()
plot_ont_lvl <- function(results = load_example_results(),
                         p2g = HPOExplorer::load_phenotype_to_genes(),
                         x_vars = c("genes",
                                    "celltypes",
                                    "log(abs(fold_change))")){

  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  requireNamespace("ggpubr")

  gene_symbol <- celltypes <- CellType <- ontLvl <- NULL;

  results[,celltypes:=length(unique(CellType[q<0.05])),by="hpo_id"]
  pcount <- p2g[,list(genes=length(unique(gene_symbol))),
                by="hpo_id"]
  pcount <- HPOExplorer::add_ont_lvl(pcount)
  r2 <- merge(results[,c("hpo_id","CellType","celltypes",
                         "p","q","fold_change")] |>
                unique(),
              pcount,
              by="hpo_id")

  plt <- function(x_var="log(genes)",
                  y_var="log(abs(fold_change))",
                  geom="hex",
                  method = "loess",
                  direction = 1,
                  ...){
    r2[,mean:=mean(get(gsub("[(]|[)]|log|abs","",y_var))),by="hpo_id"]
    gp <- ggplot(r2,aes_string(x=x_var,y=y_var)) +
      scale_fill_viridis_c(option = "plasma",
                           direction = direction)
    if(geom=="hex"){
      gp <- gp + geom_hex(...)
    } else if(geom=="violin"){
      gp <- gp + geom_violin(orientation = "y",
                             aes(fill=ontLvl),
                             ...)
    } else if(geom=="boxplot"){
      gp <- gp + geom_boxplot(orientation = "y",
                              aes(group=ontLvl,
                                  fill=mean),
                              ...)
    }else {
      gp <- gp + geom_jitter(width = 0,
                             alpha=.25,
                            ...)
    }
    gp <- gp +
      geom_smooth(method = method) +
      ggpubr::stat_cor(method = "pearson",
                       label.x.npc = .5,
                       label.y.npc = 1) +
      theme_bw()
    return(gp)
  }
  # plts1 <- lapply(c("log(genes)","log(celltypes)"), plt) |>
  #   patchwork::wrap_plots(ncol = 1)
  plts2 <- lapply(x_vars, plt,
                  y="ontLvl",
                  geom="boxplot",
                  direction = -1,
                  notch=TRUE) |>
    patchwork::wrap_plots(ncol = 1)
  return(list(data=r2,
              plot=plts2))
}


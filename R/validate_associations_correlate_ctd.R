#' Validate associations: correlate CellTypeDataset
#'
#' Validate phenotype-cell type associations can making comparisons within and
#' between results from different CellTypeDatasets.
#' @param group_var A character string specifying the column in \code{results}
#' that contains the group variable to compare across.
#' @param celltype_var A character string specifying the column in \code{results}
#' that contains the cell type variable to compare across.
#' @param hpo_agg_lvl Aggregate the data to a specific HPO ancestor level
#' during plotting to reduce figure size.
#' @param downsample Downsample the data to this many points when plotting.
#' @param add_density Add 2D density plots to the scatter plots.
#' @param rasterize_points Whether to rasterize the points in the scatter plots.
#' @param dpi Resolution of image after rasterization.
#' @param ... Additional arguments passed to \code{plot_density_cor}.
#' @inheritParams prioritise_targets
#' @inheritParams filter_ggstatsplot_subtitle
#' @inheritParams ggstatsplot::ggscatterstats
#' @inheritParams KGExplorer::filter_dt
#' @inheritParams ggrastr::rasterize
#' @export
#' @import data.table
#' @importFrom stats as.formula complete.cases
#' @examples
#' results <- load_example_results()[,.SD[seq(10000)],by="ctd"]
#' #### Across CTD ####
#' out1 <- validate_associations_correlate_ctd(results=results,
#'                                             group_var="ctd")
#' #### Within CTD: across developmental stages ####
#' filters <- list(ctd=c("HumanCellLandscape"), stage=c("Fetus","Adult"))
#' out2 <- validate_associations_correlate_ctd(results=results,
#'                                             filters=filters,
#'                                             group_var="stage")
validate_associations_correlate_ctd <- function(results=load_example_results(),
                                                filters=NULL,
                                                hpo_agg_lvl=3,
                                                hpo=HPOExplorer::get_hpo(),
                                                group_var="ctd",
                                                celltype_var="cl_name",
                                                q_threshold=0.05,
                                                downsample=NULL,
                                                add_density=TRUE,
                                                marginal=FALSE,
                                                rasterize_points=TRUE,
                                                dpi=100,
                                                point.args = list(alpha=.5),
                                                stats_idx=c(1,3,4,7),
                                                ...){
  if(rasterize_points) requireNamespace("ggrastr")
  test_id <- NULL;
  results <- map_celltype(results)
  results <- add_logfc(results)
  results <- KGExplorer::filter_dt(results,
                                   filters = filters)
  group_values <- unique(results[[group_var]])
  value.var <- intersect(c("p","q","logFC","estimate"),
                         names(results))


  messager("Casting results.")
  res2 <- results |>
    data.table::dcast.data.table(
      formula = stats::as.formula(paste0("hpo_id+",celltype_var,"~",group_var)),
      fun.aggregate = mean,
      drop = TRUE,
      value.var = value.var)
  res2 <- res2[stats::complete.cases(res2)][,test_id:=.I]

  ### All results ####
  n_celltypes.all <- length(unique(res2$cl_name))
  n_phenotypes.all <- length(unique(res2$hpo_id))
  messager(n_celltypes.all,"comparable celltypes.")
  messager(n_phenotypes.all,"comparable phenotypes.")

  #### Significant results in both groups ####
  res_sig <- res2[get(paste0("q_",group_values[1]))<q_threshold &
                  get(paste0("q_",group_values[2]))<q_threshold]

  n_celltypes.significant <- length(unique(res_sig$cl_name))
  n_phenotypes.significant <- length(unique(res_sig$hpo_id))
  messager(n_celltypes.significant,"comparable celltypes",
           paste0("@FDR<",q_threshold,"."))
  messager(n_phenotypes.significant,"comparable phenotypes",
           paste0("@FDR<",q_threshold,"."))

  #### Plot ####
  messager("Generating plots.")


  if(!is.null(hpo_agg_lvl)){
    res2 <- HPOExplorer::add_ancestor(res2,
                                      lvl = hpo_agg_lvl,
                                      hpo=hpo,
                                      force_new=TRUE)
    res_sig <- HPOExplorer::add_ancestor(res_sig,
                                         lvl = hpo_agg_lvl,
                                         hpo=hpo,
                                         force_new=TRUE)
    agg_var <- "ancestor_name"
  } else{
    agg_var <- NULL
  }


  plots <- list()
  plots[["p.all"]] <- plot_density_cor(res2,
                                       x=paste0("p_",group_values[1]),
                                       y=paste0("p_",group_values[2]),
                                       agg_var=agg_var,
                                       downsample=downsample,
                                       point.args=point.args,
                                       stats_idx=stats_idx,
                                       add_density=add_density,
                                       rasterize_points=rasterize_points,
                                       dpi=dpi,
                                       marginal=marginal
                                       )
  plots[["logFC.all"]] <- plot_density_cor(res2,
                                           x=paste0("logFC_",group_values[1]),
                                           y=paste0("logFC_",group_values[2]),
                                           downsample=downsample,
                                           agg_var=agg_var,
                                           point.args=point.args,
                                           stats_idx=stats_idx,
                                           add_density=add_density,
                                           rasterize_points=rasterize_points,
                                           dpi=dpi,
                                           marginal=marginal
                                           )
  plots[["estimate.all"]] <- plot_density_cor(res2,
                                           x=paste0("estimate_",group_values[1]),
                                           y=paste0("estimate_",group_values[2]),
                                           downsample=downsample,
                                           agg_var=agg_var,
                                           point.args=point.args,
                                           stats_idx=stats_idx,
                                           add_density=add_density,
                                           rasterize_points=rasterize_points,
                                           dpi=dpi,
                                           marginal=marginal
                                           )
  plots[["p.significant"]] <- plot_density_cor(res_sig,
                                               x=paste0("p_",group_values[1]),
                                               y=paste0("p_",group_values[2]),
                                               downsample=downsample,
                                               agg_var=agg_var,
                                               point.args=point.args,
                                               stats_idx=stats_idx,
                                               add_density=add_density,
                                               rasterize_points=rasterize_points,
                                               dpi=dpi,
                                               marginal=marginal
                                               )
  plots[["logFC.significant"]] <- plot_density_cor(res_sig,
                                                   x=paste0("logFC_",group_values[1]),
                                                   y=paste0("logFC_",group_values[2]),
                                                   downsample=downsample,
                                                   agg_var=agg_var,
                                                   point.args=point.args,
                                                   stats_idx=stats_idx,
                                                   add_density=add_density,
                                                   rasterize_points=rasterize_points,
                                                   dpi=dpi,
                                                   marginal=marginal
                                                   )
  plots[["estimate.significant"]] <- plot_density_cor(res_sig,
                                                   x=paste0("estimate_",group_values[1]),
                                                   y=paste0("estimate_",group_values[2]),
                                                   downsample=downsample,
                                                   agg_var=agg_var,
                                                   point.args=point.args,
                                                   stats_idx=stats_idx,
                                                   add_density=add_density,
                                                   rasterize_points=rasterize_points,
                                                   dpi=dpi,
                                                   marginal=marginal
  )
  messager("Gathering statistics.")
  data_stats <- lapply(plots, get_ggstatsplot_stats)
  #### Return ####
  return(list(
    plot=plots,
    data=list(all=res2,
              significant=res_sig),
    data_stats=data_stats,
    counts=list(
      n_celltypes.all=n_celltypes.all,
      n_phenotypes.all=n_phenotypes.all,
      n_celltypes.significant=n_celltypes.significant,
      n_phenotypes.significant=n_phenotypes.significant
    )
  ))
}

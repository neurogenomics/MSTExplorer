#' Validate associations: correlate CellTypeDataset
#'
#' Validate phenotype-cell type associations can making comparisons within and
#' between results from different CellTypeDatasets.
#' @param group_var A character string specifying the column in \code{results}
#' that contains the group variable to compare across.
#' @param celltype_var A character string specifying the column in \code{results}
#' that contains the cell type variable to compare across.
#' @param
#' @inheritParams prioritise_targets
#' @inheritParams KGExplorer::filter_dt
#' @export
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
                                                group_var="ctd",
                                                celltype_var="cl_name",
                                                q_threshold=0.05){
  test_id <- NULL;
  results <- map_celltype(results)
  add_logfc(results)
  results <- KGExplorer::filter_dt(results,
                                   filters = filters)
  group_values <- unique(results[[group_var]])
  value.var <- intersect(c("p","q","logFC","estimate"),
                         names(results))
  messager("Casting results.")
  res2 <- results |>
    data.table::dcast.data.table(
      formula = as.formula(paste0("hpo_id+",celltype_var,"~",group_var)),
      fun.aggregate = mean,
      drop = TRUE,
      value.var = value.var)
  res2 <- res2[complete.cases(res2)][,test_id:=.I]
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
  plots <- list()
  plots[["p.all"]] <- plot_density_cor(res2,
                                       x=paste0("p_",group_values[1]),
                                       y=paste0("p_",group_values[2])
                                       )
  plots[["logFC.all"]] <- plot_density_cor(res2,
                                           x=paste0("logFC_",group_values[1]),
                                           y=paste0("logFC_",group_values[2])
                                           )
  plots[["p.significant"]] <- plot_density_cor(res_sig,
                                               x=paste0("p_",group_values[1]),
                                               y=paste0("p_",group_values[2]))
  plots[["logFC.significant"]] <- plot_density_cor(res_sig,
                                                   x=paste0("logFC_",group_values[1]),
                                                   y=paste0("logFC_",group_values[2])
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

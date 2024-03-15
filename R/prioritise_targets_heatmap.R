prioritise_targets_heatmap <- function(top_targets,
                                       group_vars = c("CellType",
                                                      "hpo_name"),
                                       value_var = "effect"){

  requireNamespace("ComplexHeatmap")

  X_df <- top_targets |>
    data.table::dcast.data.table(formula = paste(
      paste(group_vars,collapse = " + "),"~ gene_symbol"
    ),
    fun.aggregate = mean,
    value.var = value_var,
    fill = 0,
    na.rm=TRUE)
  data.table::setorderv(X_df,cols = group_vars)
  X <- as.matrix(X_df[,-group_vars, with=FALSE])
  df <- X_df[,group_vars, with=FALSE]

  # hm <- heatmaply::heatmaply(X_df)
  n_vals <- lapply(df, function(x){length(unique(x))})
  # n_cols <- lapply(pals:::syspals,length)
  # valid_cols <- sort(unlist(n_cols[n_cols>max(unlist(n_vals))]))
  ra <- ComplexHeatmap::HeatmapAnnotation(df = df,
                          which = "row",
                          col = list(
                            CellType=stats::setNames(
                              pals::tol.rainbow(n_vals$CellType),
                              unique(df$CellType)),
                            Phenotype=stats::setNames(
                              pals::stepped2(n_vals$hpo_name),
                              unique(df$hpo_name))
                          )
  )
  ht_list <- ComplexHeatmap::Heatmap(matrix = X,
                     name = stringr::str_to_sentence(gsub("_"," ",value_var)),
                     cluster_rows = FALSE,
                     col = pals::inferno(n = 100),
                     right_annotation = ra,
                     row_dend_side = "right"
  )
  return(ht_list)
}

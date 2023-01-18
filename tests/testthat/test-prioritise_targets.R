test_that("prioritise_targets works", {

  results <- load_example_results()
  ctd <- load_example_ctd()

  #### Top only ####
  top_targets <- prioritise_targets(results = results,
                                    ctd = ctd)
  testthat::expect_equal(nrow(top_targets), 96)

  #### All results ####
  df <- prioritise_targets(results = results,
                           ctd = ctd,
                           top_n = NULL)
  testthat::expect_equal(nrow(df), 453)

  # df_agg <- dplyr::group_by(df, HPO_ID,Phenotype) |>
  #   dplyr::summarise(n_genes=data.table::uniqueN(Gene),
  #                    genes=list(unique(Gene)),
  #                    n_celltypes=length(unique(CellType)),
  #                    celltypes=list(unique(CellType))
  #                    ) |>
  #   dplyr::arrange(n_genes, celltypes) |>
  #   data.table::data.table()
  #  data.table::fwrite(df_agg,"~/Downloads/df_agg.csv")
  #
  # df_intel <- df[disease_characteristic=="Intellectual disability" & (!Phenotype %in% c("Choreoathetosis","Coma")),]
  # top_genes <- sort(table(df_intel$Gene),decreasing = TRUE)
#
#     sort(table(unique(df_intel[,c("Phenotype","HPO_ID","CellType")])$CellType),
#          decreasing = TRUE)

})

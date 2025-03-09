test_that("plot_differential_outcomes works", {

  run_tests <- function(p){
    testthat::expect_true(methods::is(p$plot[[1]],"gg"))
    testthat::expect_true(methods::is(p$data$summary_data,"data.table"))
    testthat::expect_true(methods::is(p$data$pairwise_data,"data.table"))
  }

  #### Load example data ####
  results <- load_example_results()

  results <- add_symptom_results(results = results,
                                 celltype_col="cl_name",
                                 proportion_driver_genes_symptom_threshold=.25)

  #### Multiple phenotypes per disease #####
  results <- HPOExplorer::add_gpt_annotations(results)
  annot <- HPOExplorer::load_phenotype_to_genes(3)

  # Get hypotonia-related subset of data
  diseases_with_hypotonia <- annot[hpo_id=="HP:0001319"]$disease_id |>unique()## Neonatal hypotonia
  phenotypes_with_hypotonia <- annot[disease_id %in% diseases_with_hypotonia]$hpo_id |> unique()
  results[,group:="Hypotonia"]

  # Limit the number of cell types
  # (or else there's an insane number of pairwise comparisons)
  ct_select <- unique(results$CellType)[seq(5)]

  p1a <- plot_differential_outcomes(results[q<0.05 & CellType %in% ct_select],
                                    # max_facets = 10,
                                    filters = list(disease_id=diseases_with_hypotonia),
                                    facet_var = "group",
                                    y_var = "death",
                                    run_stats = TRUE)
  sig_res <- p1a$data$summary_data[q.value<0.05]$facet
  plt <- p1a$plot[[sig_res[1]]]

  run_tests(p1a)

  p1b <- plot_differential_outcomes(results[CellType %in% ct_select],
                                    max_facets=10,
                                   facet_var = "disease_name",
                                   y_var = "blindness",
                                   run_stats = TRUE)
  run_tests(p1b)

  p1c <- plot_differential_outcomes(results[CellType %in% ct_select],
                                    max_facets=10,
                                    facet_var = "disease_name",
                                    y_var = "cancer",
                                    run_stats = TRUE)
  patchwork::wrap_plots(p1c$plot[1])
  run_tests(p1c)

  #### Multiple diseases per phenotype: Severity_score #####
  results <- HPOExplorer::add_severity(results)
  # Testthat interprets messages with errors as failures even with tryCatch()
  testthat::expect_error(
    p2 <- plot_differential_outcomes(results,
                                     max_facets=10,
                                     facet_var = "hpo_name",
                                     y_var = "Severity_score",
                                     run_stats = TRUE)
  )

  #### Multiple diseases per phenotype: AgeOfDeath #####
  results <- HPOExplorer::add_death(results,
                                    agg_by = c("disease_id","hpo_id"))
  p3 <- plot_differential_outcomes(results,
                                   filters = list(hpo_name=c(
                                     "Hypotonia",
                                     "Neonatal hypotonia"
                                   )),
                                   facet_var = "hpo_name",
                                   # x_var = "CellType",
                                   y_var = "AgeOfDeath_score_mean",
                                   run_stats = TRUE)

  run_tests(p3)


  #### Multiple diseases per phenotype: onset_score_mean #####
  results <- HPOExplorer::add_onset(results,
                                    agg_by = c("disease_id","hpo_id"))
  testthat::expect_warning(
    p4 <- plot_differential_outcomes(results,
                                     filters = list(hpo_name=c(
                                       "Hypotonia",
                                       "Neonatal hypotonia"
                                     )),
                                     facet_var = "hpo_name",
                                     y_var = "onset_score_mean",
                                     run_stats = TRUE)
  )

  run_tests(p4)
})

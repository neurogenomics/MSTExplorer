test_that("plot_differential_outcomes works", {

  run_tests <- function(p){
    testthat::expect_true(methods::is(p$plot[[1]],"gg"))
    testthat::expect_true(methods::is(p$data$summary_data,"data.table"))
    testthat::expect_true(methods::is(p$data$pairwise_data,"data.table"))
  }
  # results <- load_example_results("phenomix_results.tsv.gz")
  results <- load_example_results()
  results <- add_symptom_results(results = results,
                                 celltype_col="cl_name",
                                 proportion_driver_genes_symptom_threshold=.25)
  #### Multiple phenotypes per disease #####
  results <- HPOExplorer::add_gpt_annotations(results)
  p1a <- plot_differential_outcomes(results,
                                    max_facets=10,
                                     facet_var = "disease_name",
                                     y_var = "severity_score_gpt",
                                     run_stats = TRUE)
  run_tests(p1a)

  p1b <- plot_differential_outcomes(results,
                                    max_facets=10,
                                   facet_var = "disease_name",
                                   y_var = "blindness",
                                   run_stats = TRUE)
  run_tests(p1b)

  p1c <- plot_differential_outcomes(results,
                                    max_facets=10,
                                    facet_var = "disease_name",
                                    y_var = "cancer",
                                    run_stats = TRUE)
  patchwork::wrap_plots(p1c$plot$`Granulomatosis with polyangiitis`)
  run_tests(p1c)

  # #### Multiple diseases per phenotype: Severity_score #####
  results <- HPOExplorer::add_severity(results)
  p2 <- plot_differential_outcomes(results,
                                   max_facets=10,
                                   facet_var = "hpo_name",
                                   y_var = "Severity_score",
                                   run_stats = TRUE)

  # #### Multiple diseases per phenotype: AgeOfDeath #####
  results <- HPOExplorer::add_death(results,
                                    agg_by = c("disease_id","hpo_id"))
  p3 <- plot_differential_outcomes(results,
                                   filters = list(hpo_name=c(
                                     "Hypotonia",
                                     "Neonatal hypotonia"
                                   )),
                                   facet_var = "hpo_name",
                                   y_var = "AgeOfDeath_score_mean",
                                   run_stats = TRUE)

  # #### Multiple diseases per phenotype: onset_score_mean #####
  results <- HPOExplorer::add_onset(results,
                                    agg_by = c("disease_id","hpo_id"))
  p4 <- plot_differential_outcomes(results,
                                   filters = list(hpo_name=c(
                                     "Hypotonia",
                                     "Neonatal hypotonia"
                                   )),
                                   facet_var = "hpo_name",
                                   y_var = "onset_score_mean",
                                   run_stats = TRUE)
})

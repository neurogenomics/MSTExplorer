test_that("plot_differential_outcomes works", {

  run_tests <- function(p){
    testthat::expect_true(methods::is(p$plot[[1]],"gg"))
    testthat::expect_true(methods::is(p$data$summary_data,"data.table"))
    testthat::expect_true(methods::is(p$data$pairwise_data,"data.table"))
  }
  # results <- load_example_results("phenomix_results.tsv.gz")
  results <- load_example_results()
  # results <- map_celltype(results)
  results <- add_symptom_results(results = results,
                                 celltype_col="cl_name",
                                 proportion_driver_genes_symptom_threshold=.25)
  #### Multiple phenotypes per disease #####
  results <- HPOExplorer::add_gpt_annotations(
    results
    # annot = HPOExplorer::gpt_annot_codify(reset_tiers_dict=TRUE)$annot
    )
  annot <- HPOExplorer::load_phenotype_to_genes(3)
  diseases_with_hypotonia <- annot[hpo_id=="HP:0001319"]$disease_id |>unique()## Neonatal hypotonia
  phenotypes_with_hypotonia <- annot[disease_id %in% diseases_with_hypotonia]$hpo_id |> unique()
  # hypotonias <- grep("hypotonia",
  #                    unique(results$disease_name), ignore.case = TRUE, value = TRUE)|>
  #   grep(pattern="infantile|neonatal|congenital",
  #        ignore.case = TRUE, value = TRUE)
  results[,group:="Hypotonia"]
  p1a <- plot_differential_outcomes(results[q<0.05 ],
                                    # max_facets = 10,
                                    filters = list(disease_id=diseases_with_hypotonia),
                                    facet_var = "group",
                                    y_var = "death",
                                    run_stats = TRUE)
  sig_res <- p1a$data$summary_data[q.value<0.05]$facet
  plt <- p1a$plot[[sig_res[1]]] #+
    ggplot2::scale_y_continuous(labels = c("never (0)","rarely (1)","variable (2)","often (3)"))
  # ggplot2::ggsave(filename = "~/Downloads/devdelay_hypotonia.png",plt, height=5, width=7)
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
                                   # x_var = "CellType",
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

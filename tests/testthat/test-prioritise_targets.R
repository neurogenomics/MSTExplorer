test_that("prioritise_targets works", {

  results <- load_example_results()[seq(50000),]
  ctd_list <- load_example_ctd(c("ctd_DescartesHuman.rds",
                                 "ctd_HumanCellLandscape.rds"),
                                multi_dataset = TRUE)

  #### Top only ####
  res1 <- prioritise_targets(results = results,
                             ctd_list = ctd_list)
  testthat::expect_gte(nrow(res1$top_targets), 6)

  #### All results ####
  res2 <- prioritise_targets(results = results,
                             ctd_list = ctd_list,
                             top_n = 2)
  testthat::expect_gte(nrow(res2$top_targets), 6)

  #### Plot evidence score vs. specificity ####
  res3 <- prioritise_targets(results = results[seq(2500)],
                             ctd_list = ctd_list,
                             keep_deaths = NULL,
                             #### Phenotype level ####
                             keep_ont_levels = NULL,
                             pheno_ndiseases_threshold = NULL,
                             keep_tiers = NULL,
                             severity_threshold_max = NULL,
                             #### Symptom level ####
                             severity_threshold = NULL,
                             pheno_frequency_threshold = NULL,
                             keep_onsets = NULL,
                             #### Celltype level ####
                             q_threshold = 0.05,
                             effect_threshold = 1,
                             symptom_p_threshold = NULL,
                             symptom_intersection_size_threshold = 1,
                             keep_celltypes = NULL,
                             #### Gene level ####
                             evidence_score_threshold = NULL,
                             keep_chr = NULL,
                             gene_size = list("min"=0,
                                              "max"=Inf),
                             gene_frequency_threshold = NULL,
                             keep_biotypes = NULL,
                             keep_specificity_quantiles = NULL,
                             keep_mean_exp_quantiles = seq(1,40),
                             symptom_gene_overlap = TRUE)
  testthat::expect_gte(nrow(res3$top_targets), 40000)
})

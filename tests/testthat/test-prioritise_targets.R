test_that("prioritise_targets works", {

  testthat::skip() # skip for now

  # Tag the specific HPO release version for reproducibility
  tag <- "2024-02-08"
  phenotype_to_genes = HPOExplorer::load_phenotype_to_genes(tag=tag)
  hpo = HPOExplorer::get_hpo(tag = tag)

  # Tag the specific MSTExplorer release version for reproducibility
  results <- load_example_results(tag = "v0.1.10" )[q<0.01]
  ctd_list <- load_example_ctd(c("ctd_DescartesHuman.rds",
                                 "ctd_HumanCellLandscape.rds"),
                                multi_dataset = TRUE)

  #### Top only ####
  res1 <- prioritise_targets(results = results,
                             phenotype_to_genes = phenotype_to_genes,
                             hpo = hpo,
                             ctd_list = ctd_list)
  testthat::expect_gte(nrow(res1$top_targets), 6)

  #### All results ####
  res2 <- prioritise_targets(results = results,
                             phenotype_to_genes = phenotype_to_genes,
                             hpo = hpo,
                             ctd_list = ctd_list,
                             top_n = 2)
  testthat::expect_gte(nrow(res2$top_targets), 6)

  #### Plot evidence score vs. specificity ####
  res3 <- prioritise_targets(results = results,
                             phenotype_to_genes = phenotype_to_genes,
                             hpo = hpo,
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
                             keep_celltypes = NULL,
                             #### Gene level ####
                             evidence_score_threshold = NULL,
                             keep_chr = NULL,
                             gene_size = list("min"=0,
                                              "max"=Inf),
                             gene_frequency_threshold = NULL,
                             keep_biotypes = NULL,
                             keep_specificity_quantiles = NULL,
                             keep_mean_exp_quantiles = seq(1,40))
  testthat::expect_gte(nrow(res3$top_targets), 1000)
})

test_that("gen_results works", {

  set.seed(2023)

  gene_data <- HPOExplorer::load_phenotype_to_genes()
  ctd <- load_example_ctd()
  list_names <- unique(gene_data$Phenotype)[seq_len(3)]
  all_results <- gen_results(ctd = ctd,
                             gene_data = gene_data,
                             list_names = list_names,
                             reps = 10)
  testthat::expect_true(methods::is(all_results,"data.table"))
  testthat::expect_true(all(list_names %in% unique(all_results$Phenotype)))
  testthat::expect_gte(nrow(all_results[q<=0.05,]),40)
})

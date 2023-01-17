test_that("count_results works", {

  results <- load_example_results()
  celltype_counts <- count_results(results=results, group_var="CellType")
  phenotype_counts <- count_results(results=results, group_var="Phenotype")

  testthat::expect_lte(nrow(celltype_counts),
                         length(unique(results$CellType)))
  testthat::expect_lte(nrow(phenotype_counts),
                       length(unique(results$Phenotype)))
})

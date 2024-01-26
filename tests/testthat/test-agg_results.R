test_that("agg_results works", {

  phenos <- subset_results(filters=list(CellType = "Microglia"))
  agg_res <- agg_results(phenos = phenos)
  testthat::expect_true(nrow(agg_res)==1)
  testthat::expect_true(nrow(agg_res)<nrow(phenos))
})

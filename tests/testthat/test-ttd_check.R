test_that("ttd_check works", {

  top_targets <- MultiEWCE::example_targets$top_targets
  top_targets <- HPOExplorer::add_genes(top_targets[q<0.05])
  res <- ttd_check(top_targets=top_targets)
  testthat::expect_gte(nrow(res$data),7000)
  testthat::expect_gte(nrow(res$data_overlap),3000)
  testthat::expect_true(methods::is(res$plot,"ggplot"))
})

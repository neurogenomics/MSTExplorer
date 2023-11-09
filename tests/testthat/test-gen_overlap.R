test_that("gen_overlap works", {

  gene_data <- HPOExplorer::load_phenotype_to_genes()
  list_names <- unique(gene_data$disease_id)[seq(3)]
  overlap <- gen_overlap(gene_data = gene_data,
                         list_names = list_names)
  testthat::expect_gte(nrow(overlap),200)
})

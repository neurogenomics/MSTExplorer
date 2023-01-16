test_that("multiplication works", {

  gene_data <- HPOExplorer::load_phenotype_to_genes()
  gene_list <- get_gene_list(list_name="Focal non-motor seizure",
                             gene_data = gene_data)
  testthat::expect_length(gene_list,69)
})

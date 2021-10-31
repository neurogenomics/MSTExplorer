skip("dont run")
library(HPOExplorer)

if (!dir.exists("data")) {
  dir.create("data")
}


phenotype_to_genes = load_phenotype_to_genes("data/phenotype_to_genes.txt")
load("data/CTD_DescartesHuman.rda")
resultsdir = "data/results"

if (!dir.exists(resultsdir)) {
  dir.create(resultsdir)
}

Phenos1 = unique(phenotype_to_genes$Phenotype)[1:5]
results1 <- ewce_para(Phenos1,
                      phenotype_to_genes,
                      list_name_column = "Phenotype",
                      gene_column ="Gene",
                      results_directory = resultsdir,
                      ctd_file = ctd,
                      background_genes = unique(phenotype_to_genes$Gene),
                      bootstrap_reps = 10,
                      annotation_Level = 1,
                      genes_Species = "human",
                      ctd_Species = "human",
                      cores = 1)


test_that("First batch of 5 results were generated", {
  expect_equal(length(list.files(resultsdir)), 5)
  expect_equal(length(unique(merge_results(resultsdir, list_name_column = "Phenotype")$Phenotype)), 5)
})

for (f in list.files(paste0(resultsdir))) {
  file.remove(paste0(resultsdir,"/",f))
}

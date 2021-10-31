skip("Dont run")
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
Phenos2 = unique(phenotype_to_genes$Phenotype)[1:7]

results_merged <- gen_results(
  ctd,
  gene_data=phenotype_to_genes,
  list_names = Phenos1,
  background_genes = unique(phenotype_to_genes$Gene),
  list_name_column = "Phenotype",
  gene_column = "Gene",
  results_dir = resultsdir,
  overwrite_past_analysis = FALSE,
  reps = 10,
  annotLevel = 1,
  genelistSpecies = "human",
  sctSpecies = "human",
  cores = 1,
  MergeResults = FALSE
)

test_that("gen_results returns merged results and individual results in results directory",{
  expect_equal(length(list.files(resultsdir)), 5)
  expect_equal(length(unique(results_merged$Phenotype)), 5)
})

results_merged <- gen_results(
  ctd,
  gene_data=phenotype_to_genes,
  list_names = Phenos2,
  background_genes = unique(phenotype_to_genes$Gene),
  list_name_column = "Phenotype",
  gene_column = "Gene",
  results_dir = resultsdir,
  overwrite_past_analysis = FALSE,
  reps = 10,
  annotLevel = 1,
  genelistSpecies = "human",
  sctSpecies = "human",
  cores = 1,
  MergeResults = FALSE
)

test_that("gen_results can resume analysis after a stop", {
  expect_equal(length(list.files(resultsdir)), 7)
  expect_equal(length(unique(results_merged$Phenotype)), 7)
  expect_equal(length(unique(results_merged$Phenotype)), 7)
})

# Remove generated results from test
for (f in list.files(paste0(resultsdir))) {
  file.remove(paste0(resultsdir,"/",f))
}

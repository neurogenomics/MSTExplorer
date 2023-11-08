test_that("ewce_para works", {

  set.seed(2023)

  gene_data <- HPOExplorer::load_phenotype_to_genes()
  ctd <- MultiEWCE::load_example_ctd()
  list_names <- unique(gene_data$hpo_id)[seq(5)]

  #### Return results directly ####
  res_files <- ewce_para(ctd = ctd,
                         gene_data = gene_data,
                         list_names = list_names,
                         reps = 10,
                         save_dir_tmp = NULL)
  for(r in res_files){
    testthat::expect_true(
      all(c("results","hit.cells","bootstrap_data") %in% names(r))
    )
  }
  #### Return paths to saved results ####
  save_dir_tmp <- file.path(tempdir(),"results")
  res_files2 <- ewce_para(ctd = ctd,
                          gene_data = gene_data,
                          list_names = list_names,
                          reps = 10,
                          save_dir_tmp = save_dir_tmp)
  files <- list.files(save_dir_tmp, full.names = TRUE)
  testthat::expect_lte(length(files), length(list_names))
  for (f in res_files2) {
   testthat::expect_true(file.exists(f))
   r2 <- readRDS(f)
   testthat::expect_true(all(c("results","hit.cells","bootstrap_data") %in% names(r2)))
  }

  #### Tests get_unfinished_list_names ####
  all_phenotypes <- unique(gene_data$hpo_id)
  unfinished <- get_unfinished_list_names(list_names = all_phenotypes,
                                          save_dir_tmp = save_dir_tmp)
  testthat::expect_lte(length(unfinished),
                       length(all_phenotypes))


  #### Merge results ####
  all_results1 <- merge_results(res_files=res_files)
  all_results2 <- merge_results(res_files=res_files2)
  ## Confirm both methods have the correct phenotyoes
  testthat::expect_gte(sum(list_names %in% unique(all_results1$hpo_id)),
                       length(list_names)-3)
  testthat::expect_gte(sum(list_names %in% unique(all_results2$hpo_id)),
                       length(list_names)-3)
  ## Confirm both methods are identical
  data.table::setkey(all_results1,"CellType")
  data.table::setkey(all_results2,"CellType")
  testthat::expect_equal(nrow(all_results1),
                         nrow(all_results2))
  testthat::expect_equal(sort(unique(all_results1$hpo_id)),
                         sort(unique(all_results2$hpo_id)))
  testthat::expect_equal(sort(unique(all_results1$CellType)),
                         sort(unique(all_results2$CellType)))
  # testthat::expect_equal(all_results1, all_results2)
})

#!/usr/bin/env Rscript
library("optparse")

option_list <-list(
  optparse::make_option(c("-i", "--idx"), type="integer", default=1,
                        help="PBS_ARRAY_INDEX", metavar="character"),
  optparse::make_option( c("-n", "--ncpus"), type="integer", default=4,
                         help="Number of CPUs to use.", metavar="character"),
  optparse::make_option( c("-b", "--batches"), type="integer", default=100,
                         help="Number of total batches.", metavar="character")
);
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)
root <- "/rds/general/project/neurogenomics-lab/ephemeral/rare_disease"


library(MSTExplorer)
ctd <- load_example_ctd("ctd_DescartesHuman.rds")
gene_data <- HPOExplorer::load_phenotype_to_genes()
gene_data[,n_gene:=length(unique(gene_symbol)),by="hpo_id"]
gene_data <- gene_data[n_gene>=4,]
#### Split HPO IDs into N chunks ####
ids <- unique(gene_data$hpo_id)
chunks <- split(ids, cut(seq_along(ids),opt$batches,labels = FALSE))

#### Run enrichment analyses #####
all_results <- gen_results(ctd = ctd,
                           list_name_column = "hpo_id",
                           list_names = chunks[[opt$idx]],
                           gene_data = gene_data,
                           annotLevel = 2,
                           reps = 100000,
                           cores = opt$ncpus,
                           save_dir = file.path(root,paste0("batch",opt$idx)))


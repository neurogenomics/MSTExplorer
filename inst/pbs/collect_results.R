
# root <- "/rds/general/project/neurogenomics-lab/ephemeral/rare_disease"
root <- "/Volumes/bms20/projects/neurogenomics-lab/ephemeral/rare_disease"
f <- list.files(root, full.names = TRUE, pattern = ".rds", recursive = TRUE)
res <- lapply(f,readRDS) |>
  `names<-`(basename(dirname(f))) |>
  data.table::rbindlist(use.names = TRUE,idcol = "batch", fill = TRUE)
length(unique(res$hpo_id))
length(unique(res$hpo_id))*length(unique(res$CellType))
#### Ensure only one test per celltype-phenotype combination ####
res <- res[,utils::head(.SD, 1),by = c("hpo_id","CellType")]
## Recalculate q across all batches
res$q <- stats::p.adjust(res$p,method = "bonf")
res[q<0.05,]
length(unique(res[q<0.05,]$hpo_id))

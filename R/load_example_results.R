#' Load example results
#'
#' This loads a example of enrichment results from \link[MultiEWCE]{gen_results}
#' using the Human Phenotype Ontology and a given CellTypeDataset (CTD).
#' Results were then merged together with  \link[MultiEWCE]{merge_results}.
#' @param file File to load:
#' \itemize{
#' \item{"Descartes_All_Results_extras.symptoms.full_join.rds":}{
#' Contains cell type-phenotype and cell type-symptom (phenotype + disease)
#' enrichment results merged into one table.
#' Used the \href{https://descartes.brotmanbaty.org/}{Descartes}
#' CTD ( annotation level 1).
#' }
#' \item{"Descartes_All_Results_extras.rds":}{
#' Contains cell type-phenotype enrichment results.
#' Used the \href{https://descartes.brotmanbaty.org/}{Descartes}
#' CTD ( annotation level 1).
#' }
#' \item{"tabulamuris_merged.rds":}{
#' Contains cell type-phenotype enrichment results.
#' Used the \href{https://tabula-muris.ds.czbiohub.org/}{Tabula Muris}
#' CTD.
#' }
#' }
#' @param force_new Download the file even when a local copy already exists.
#' @inheritParams load_example_ctd
#' @inheritParams piggyback::pb_download
#' @source
#' \code{
#' d <-  "~/Desktop/ewce/rare_disease_celltyping_apps/cell_select"
#' #### Descartes_All_Results_extras ####
#' f0 <- file.path(d,"data/Descartes_All_Results_extras.rds")
#' r0 <- readRDS(f0)
#' data.table::setnames(r0,"list","Phenotype")
#' r0$HPO_id=NULL
#' r0 <- HPOExplorer:::fix_hpo_ids(dt=r0)
#' f0new <- file.path(tempdir(),basename(f0))
#' saveRDS(r0,file = f0new)
#' piggyback::pb_upload(file = f0new,
#'                      tag = "v0.0.1", repo = "neurogenomics/MultiEWCE")
#' #### Descartes_All_Results ####
#' f1 <- file.path(d,"data/Descartes_All_Results.rds")
#' r1 <- readRDS(f1)
#' data.table::setnames(r1,"list","Phenotype")
#' r1 <- HPOExplorer:::fix_hpo_ids(dt=r1)
#' f1new <- file.path(tempdir(),basename(f1))
#' saveRDS(r1,file = f1new)
#' piggyback::pb_upload(file = f1new,
#'                      tag = "v0.0.1", repo = "neurogenomics/MultiEWCE")
#' #### tabulamuris_merged ####
#' f2 <- file.path(d,"data/tabulamuris_merged.rds")
#' r2 <- readRDS(f2)
#' data.table::setnames(r2,"list","Phenotype")
#' r2 <- HPOExplorer:::fix_hpo_ids(dt=r2)
#' f2new <- file.path(tempdir(),basename(f2))
#' saveRDS(r2,file = f2new)
#' piggyback::pb_upload(file = f2new,
#'                      tag = "v0.0.1", repo = "neurogenomics/MultiEWCE")
#' }
#' @source \href{https://github.com/neurogenomics/rare_disease_celltyping/}{
#' Results located in 'results' folder.}
#' @returns dataframe of enrichment results.
#'
#' @export
#' @importFrom piggyback pb_download
#' @importFrom tools R_user_dir
#' @examples
#' res <- load_example_results()
load_example_results <- function(file=c(
  "Descartes_All_Results_extras.symptoms.rds",
  "Descartes_All_Results_extras.symptoms.full_join.rds",
  "Descartes_All_Results_extras.rds",
  "gen_overlap.symptoms.filt.rds",
  "tabulamuris_merged.rds"),
  tag = "v0.0.1",
  save_dir=tools::R_user_dir(package = "MultiEWCE"),
  force_new=FALSE
  ) {

  file <- file[[1]]
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  save_path <- file.path(save_dir,file)
  if(file.exists(save_path) && isTRUE(force_new)){
    rm_ <- file.remove(save_path)
  }
  if (!file.exists(save_path)) {
    piggyback::pb_download(file = basename(file),
                           repo = "neurogenomics/MultiEWCE",
                           tag = tag,
                           dest = save_dir,
                           overwrite = TRUE)
  }
  results <- readRDS(save_path)
  names(results) <- gsub("^symptoms\\.","symptom.",names(results))
  return(results)
}

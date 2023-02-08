#' Load example results
#'
#' This loads a example of enrichment results from \link[MultiEWCE]{gen_results}
#' using the Human Phenotype Ontology and a given CellTypeDataset (CTD).
#' Results were then merged together with  \link[MultiEWCE]{merge_results}.
#' @param file File to load:
#' \itemize{
#' \item{"Descartes_All_Results.rds":}{
#' Used the \href{https://descartes.brotmanbaty.org/}{Descartes}
#' CTD ( annotation level 1).
#' }
#' \item{"tabulamuris_merged.rds":}{
#' Used the \href{https://tabula-muris.ds.czbiohub.org/}{Tabula Muris}
#' CTD.
#' }
#' }
#' @inheritParams load_example_ctd
#' @inheritParams piggyback::pb_download
#' @source
#' \code{
#' d <-  "~/Desktop/ewce/rare_disease_celltyping_apps/cell_select_interactive"
#' piggyback::pb_upload(file = file.path(d,"data/Descartes_All_Results.rds"),
#'                      tag = "v0.0.1", repo = "neurogenomics/MultiEWCE")
#' piggyback::pb_upload(file = file.path(d,"data/tabulamuris_merged.rds"),
#'                      tag = "v0.0.1", repo = "neurogenomics/MultiEWCE")
#' }
#' @source \href{https://github.com/neurogenomics/rare_disease_celltyping/}{
#' Results located in 'results' folder.}
#' @returns dataframe of enrichment results.
#'
#' @export
#' @importFrom piggyback pb_download
#' @importFrom tools R_user_dir
#' @importFrom data.table setkeyv
#' @examples
#' res <- load_example_results()
load_example_results <- function(file=c("Descartes_All_Results_extras.rds",
                                        "tabulamuris_merged.rds"),
                                 tag = "v0.0.1",
                                 save_dir=tools::R_user_dir(package = "MultiEWCE")
                                 ) {

  file <- file[[1]]
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  save_path <- file.path(save_dir,file)
  if (!file.exists(save_path)) {
    piggyback::pb_download(file = basename(file),
                           repo = "neurogenomics/MultiEWCE",
                           tag = tag,
                           dest = save_dir,
                           overwrite = TRUE)
  }
  results <- readRDS(save_path)
  data.table::setnames(results,"list","Phenotype")
  if(file=="tabulamuris_merged.rds"){
    results$HPO_ID <- HPOExplorer::harmonise_phenotypes(
      phenotypes = results$Phenotype,
      as_hpo_ids = TRUE)
  } else {
    results <- HPOExplorer::add_hpo_id(phenos = results,
                                       verbose = FALSE)
  }
  return(results)
}

#' Load example results
#'
#' This loads a example of enrichment results from \link[MSTExplorer]{gen_results}
#' using the Human Phenotype Ontology and a given CellTypeDataset (CTD).
#' Results were then merged together with  \link[MSTExplorer]{merge_results}.
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
#' @param multi_dataset Merge results generated from
#'  multiple CellTypeDataset references
#'  (e.g. DescartesHuman and HumanCellLandscape).
#' @param force_new Download the file even when a local copy already exists.
#' @inheritParams load_example_ctd
#' @inheritParams piggyback::pb_download
#' @source
#' \code{
#' d <-  "~/Desktop/ewce/rare_disease_celltyping_apps/cell_select"
#' #### Descartes_All_Results_extras ####
#' f0 <- file.path(d,"data/Descartes_All_Results_extras.rds")
#' r0 <- readRDS(f0)
#' data.table::setnames(r0,"list","hpo_name")
#' r0$HPO_id=NULL
#' r0 <- HPOExplorer:::fix_hpo_ids(dt=r0)
#' f0new <- file.path(tempdir(),basename(f0))
#' saveRDS(r0,file = f0new)
#' piggyback::pb_upload(file = f0new,
#'                      tag = "v0.0.1", repo = "neurogenomics/MSTExplorer")
#' #### Descartes_All_Results ####
#' f1 <- file.path(d,"data/Descartes_All_Results.rds")
#' r1 <- readRDS(f1)
#' data.table::setnames(r1,"list","hpo_name")
#' r1 <- HPOExplorer:::fix_hpo_ids(dt=r1)
#' f1new <- file.path(tempdir(),basename(f1))
#' saveRDS(r1,file = f1new)
#' piggyback::pb_upload(file = f1new,
#'                      tag = "v0.0.1", repo = "neurogenomics/MSTExplorer")
#' #### tabulamuris_merged ####
#' f2 <- file.path(d,"data/tabulamuris_merged.rds")
#' r2 <- readRDS(f2)
#' data.table::setnames(r2,"list","hpo_name")
#' r2 <- HPOExplorer:::fix_hpo_ids(dt=r2)
#' f2new <- file.path(tempdir(),basename(f2))
#' saveRDS(r2,file = f2new)
#' piggyback::pb_upload(file = f2new,
#'                      tag = "v0.0.1", repo = "neurogenomics/MSTExplorer")
#' }
#' @source \href{https://github.com/neurogenomics/rare_disease_celltyping/}{
#' Results located in 'results' folder.}
#' @returns dataframe of enrichment results.
#'
#' @export
#' @importFrom piggyback pb_download
#' @examples
#' res <- load_example_results()
load_example_results <- function(file=c(
  "phenomix_results.tsv.gz",
  # "rare_disease_min_genes4_DescartesHuman.rds",
  # "rare_disease_min_genes4_HumanCellLandscape.rds"
  # "results_DescartesHuman.csv.gz",
  # "rare_disease_min_genes4_DescartesHuman.rds",
  # "Descartes_All_Results_extras.symptoms.rds",
  # "Descartes_All_Results_extras.symptoms.full_join.rds",
  # "Descartes_All_Results_extras.rds",
  # "gen_overlap.symptoms.filt.rds",
  # "tabulamuris_merged.rds"
  ),
  multi_dataset=FALSE,
  tag = "latest",
  save_dir=KGExplorer::cache_dir(package="MSTExplorer"),
  force_new=FALSE
  ) {
  if(multi_dataset){
    ctd_names <- c("DescartesHuman","HumanCellLandscape")
    res <- lapply(stats::setNames(ctd_names,ctd_names), function(x){
      load_example_results(paste0("rare_disease_min_genes4_",x,".rds"))
    })  |> data.table::rbindlist(idcol = "ctd")
    return(res)
  } else{
    file <- file[[1]]
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    save_path <- file.path(save_dir,file)
    if(file.exists(save_path) && isTRUE(force_new)){
      rm_ <- file.remove(save_path)
    }
    if (!file.exists(save_path)) {
      piggyback::pb_download(file = basename(file),
                             repo = "neurogenomics/MSTExplorer",
                             tag = tag,
                             dest = save_dir,
                             overwrite = TRUE)
    }
    if(grepl("\\.rds$",save_path, ignore.case = TRUE)){
      results <- readRDS(save_path)
    } else {
      results <- data.table::fread(save_path)
    }
    data.table::setnames(results,
                         c("HPO_ID","Phenotype","HPO_ID.disease_id","disease_id",
                           "HPO_ID.LinkID","LinkID"),
                         c("hpo_id","hpo_name","hpo_id.disease_id","disease_id",
                           "hpo_id.disease_id","disease_id"),
                         skip_absent = TRUE)
    names(results) <- gsub("^symptoms\\.","symptom.",names(results))
    return(results)
  }

}

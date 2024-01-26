#' Map cell type
#'
#' Map cell types to cell ontology terms.
#' @param map 1:1 mapping between cell ontology terms and cell type names used
#'  in \code{results}.
#' @param rm_prefixes Prefixes to remove from cell type names in
#' \code{input_col} before performing mapping.
#' @param by Columns to merge on.
#' @inheritParams ggnetwork_plot_full
#' @export
#' @import data.table
#' @examples
#' results <- load_example_results(multi_dataset = TRUE)
#' results2 <- map_celltype(results = results)
map_celltype <- function(results,
                         input_col="CellType",
                         map = KGExplorer::get_data_package(
                           package = "MultiEWCE",
                           name="celltype_maps"),
                         rm_prefixes=c("Adult","Fetus","HESC"),
                         by=c("ctd","author_celltype")
                         ){
  author_celltype <- NULL;
  new_cols <- c("cell_type_ontology_term_id","cell_type")
  if(all(new_cols %in% names(results))) {
    return(results)
  }
  results[,author_celltype:=gsub(paste(paste0("^",rm_prefixes,"_"),
                                       collapse = "|"),"",
                                 get(input_col),ignore.case = TRUE)]
  if(!all(by %in% names(results))) {
    stopper("All 'by' columns must be in 'results'.")
  }
  results_cl <- data.table::merge.data.table(results,
                               map,
                               by=by,
                               all.x = TRUE)
  if(sum(is.na(results_cl$cell_type_ontology_term_id))>0){
    stopper("Missing 'cell_type_ontology_term_id' for",
            sum(is.na(results_cl$cell_type_ontology_term_id)),"rows.")
  }
  ## Rename cols to make more concise and conform to hpo_id/hpo_name format
  data.table::setnames(results_cl,
                       c("cell_type_ontology_term_id","cell_type"),
                       c("cl_id","cl_name"))
  #### Add stage ####
  # results[ctd=="DescartesHuman",stage:="Fetus"]
  # results[ctd=="HumanCellLandscape",
  #         stage:=data.table::tstrsplit(CellType,"_",keep=1)]
  return(results_cl)
}



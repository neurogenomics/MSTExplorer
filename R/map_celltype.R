#' Map cell type
#'
#' Map cell types to cell ontology terms.
#' @param map 1:1 mapping between cell ontology terms and cell type names used
#'  in \code{results}.
#' @param rm_prefixes Prefixes to remove from cell type names in
#' \code{input_col} before performing mapping.
#' @param by Columns to merge on.
#' @param input_col Column to use for linking with the \code{map} data.
#' @param add_stage Add developmental stage information.
#' @inheritParams ggnetwork_plot_full
#' @export
#' @import data.table
#' @examples
#' results <- load_example_results()
#' results2 <- map_celltype(results = results)
map_celltype <- function(results,
                         input_col="CellType",
                         map = KGExplorer::get_data_package(
                           package = "MSTExplorer",
                           name="celltype_maps"),
                         rm_prefixes=c("Adult","Fetus","HESC"),
                         by=c("ctd","author_celltype"),
                         add_stage=TRUE
                         ){
  author_celltype <- ctd <- CellType <- NULL;
  #### Check for existing columns ####
  new_cols <- c("cl_id","cl_name")
  if(all(new_cols %in% names(results))) {
    messager("Cell type columns already present.","Skipping mapping.")
    return(results)
  }
  #### Add new columns ####
  messager("Mapping cell types to cell ontology terms.")
  results[,author_celltype:=gsub(paste(paste0("^",rm_prefixes,"_"),
                                       collapse = "|"),"",
                                 get(input_col),ignore.case = TRUE)]
  by <- by[by %in% names(results)]
  if(length(by)==0) {
    stopper("All 'by' columns must be in 'results'.")
  }
  results_cl <- data.table::merge.data.table(results,
                               map,
                               by=by,
                               all.x = TRUE)
  if(sum(is.na(results_cl$cl_id))>0){
    stopper("Missing 'cl_id' for",
            sum(is.na(results_cl$cl_id)),"rows.")
  }
  #### Add stage ####
  if(isTRUE(add_stage)){
    messager("Adding stage information.")
    results_cl[ctd=="DescartesHuman",
               stage:="Fetus"]
    results_cl[ctd=="HumanCellLandscape",
               stage:=data.table::tstrsplit(CellType,"_",keep=1)]
  }
  return(results_cl)
}



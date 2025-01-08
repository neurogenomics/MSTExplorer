#' Map cell type
#'
#' Map cell types to cell ontology terms.
#' @param map Many:many mapping between tissue ontology terms and cell type
#' ontology terms.
#' @param return_agg Return the aggregated results instead of merging with a
#' table of association \code{results}.
#' @param uberon UBERON ontology object of class \link[simona]{ontology_DAG}.
#' @inheritParams ggnetwork_plot_full
#' @inheritParams KGExplorer::get_ontology
#' @export
#' @import data.table
#' @examples
#' results <- load_example_results()
#' results2 <- map_tissue(results = results)
map_tissue <- function(results = NULL,
                       map = KGExplorer::get_data_package(
                           package = "MSTExplorer",
                           name="tissue_maps"),
                       lvl=10,
                       uberon = KGExplorer::get_ontology(name = "uberon",
                                                         lvl=lvl),
                       return_agg=FALSE
                       ){
  ancestor <- ancestor_name <- cl_count <- uberon_id <- uberon_name <- id <-
    uberon_ancestor_name <- top_uberon_name <- uberon_ancestor <-
    top_uberon_id <- NULL;
  if(!is.null(results)){
    results <- map_celltype(results)
    new_cols <- c("top_uberon_id","uberon_ancestor_name")
    if(all(new_cols %in% names(results))) {
      messager("Tissue columns already present.","Skipping mapping.")
      return(results)
    }
  }
  messager("Mapping cell types to UBERON tissue ontology terms.")
  #### Assign each tissue to a top tissue ####
  # results[!cl_id %in% tissue_maps$cl_id]
  by <- c("ctd","cl_id","cl_name")
  if(!is.null(results)) by <- intersect(by,names(results))
  map_agg <- map[,list(
    cl_count=sum(cl_count),
    n_uberon=data.table::uniqueN(uberon_id),
    top_uberon_id=names(table(uberon_id))[which.max(table(uberon_id))],
    top_uberon_name=names(table(uberon_name))[which.max(table(uberon_name))]
                       ),
    by=by]
  #### Get the ancestor for each tissue #####
  if(!isFALSE(lvl)){
    ancestor_dat <- data.frame(uberon@elementMetadata@listData) |>
      dplyr::select(id,ancestor,ancestor_name) |>
      dplyr::rename(
        uberon_id = id,
        uberon_ancestor = ancestor,
        uberon_ancestor_name = ancestor_name
      )
    map_agg2 <- merge(map_agg,
                      ancestor_dat,
                      all.x=TRUE,
                      by.x="top_uberon_id",
                      by.y="uberon_id")
    map_agg2[,uberon_ancestor_name:=data.table::fcoalesce(uberon_ancestor_name,
                                                          top_uberon_name)]
    map_agg2[,uberon_ancestor:=data.table::fcoalesce(uberon_ancestor,
                                                     top_uberon_id)]
  } else {
    map_agg2 <- map_agg
  }
  if(isTRUE(return_agg) || is.null(results)) return(map_agg2)

  results2 <- data.table::merge.data.table(results,
                                           map_agg2,
                                           all.x = TRUE,
                                           by=by)
  return(results2)
}



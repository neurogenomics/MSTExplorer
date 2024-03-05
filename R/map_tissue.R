#' Map cell type
#'
#' Map cell types to cell ontology terms.
#' @param map Many:many mapping between tissue ontology terms and cell type
#' ontology terms.
#' @inheritParams ggnetwork_plot_full
#' @export
#' @import data.table
#' @examples
#' results <- load_example_results()
#' results2 <- map_tissue(results = results)
map_tissue <- function(results,
                         map = KGExplorer::get_data_package(
                           package = "MSTExplorer",
                           name="tissue_maps"),
                       add_ancestors=10
                       ){
  ancestor <- ancestor_name <- n_tissues <-
    n_ancestors <- cl_name <- NULL;
  results <- map_celltype(results)
  new_cols <- c("top_uberon_id","uberon_ancestor_name")
  if(all(new_cols %in% names(results))) {
    messager("Tissue columns already present.","Skipping mapping.")
    return(results)
  }
  messager("Mapping cell types to UBERON tissue ontology terms.")
  #### Assign each tissue to a top tissue ####
  # results[!cl_id %in% tissue_maps$cl_id]
  by <- c("ctd","cl_id","cl_name")
  by <- intersect(by,names(results))
  map_agg <- map[,list(cl_count=sum(cl_count),
                       n_uberon=data.table::uniqueN(uberon_id),
                       top_uberon_id=names(table(uberon_id))[which.max(table(uberon_id))],
                       top_uberon_name=names(table(uberon_name))[which.max(table(uberon_name))]
                       ),
                 by=by]
  #### Get the ancestor for each tissue #####
  uberon <- KGExplorer::get_ontology("uberon",
                                     add_ancestors=add_ancestors)
  map_agg2 <- merge(map_agg,
        data.table::data.table(
          uberon@elementMetadata
          )[,list(uberon_id=id,
                  uberon_ancestor=ancestor,
                  uberon_ancestor_name=ancestor_name)],
        all.x=TRUE,
        by.x="top_uberon_id",
        by.y="uberon_id")
  map_agg2[,uberon_ancestor_name:=data.table::fcoalesce(uberon_ancestor_name,top_uberon_name)]
  map_agg2[,uberon_ancestor:=data.table::fcoalesce(uberon_ancestor,top_uberon_id)]
  results2 <- data.table::merge.data.table(results,
                                           map_agg2,
                                           all.x = TRUE,
                                           by=by)
  return(results2)
}



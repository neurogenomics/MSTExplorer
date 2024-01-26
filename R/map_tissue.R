#' Map cell type
#'
#' Map cell types to cell ontology terms.
#' @param map Many:many mapping between tissue ontology terms and cell type
#' ontology terms.
#' @inheritParams ggnetwork_plot_full
#' @export
#' @import data.table
#' @examples
#' results <- load_example_results(multi_dataset = TRUE)
#' results2 <- map_tissue(results = results)
map_tissue <- function(results,
                         map = KGExplorer::get_data_package(
                           package = "MultiEWCE",
                           name="anatomy_maps")
                       ){
  results <- map_celltype(results)
  new_cols <- c("tissue_ontology_term_id","tissue")
  if(all(new_cols %in% names(results))) {
    return(results)
  }
  cl <- KGExplorer::get_ontology("cl")
  #### Assign each tissue an ancestral group ####
  uberon <- KGExplorer::get_ontology("uberon",
                                     term = map$tissue_ontology_term_id,
                                     add_ancestors = 1)
  map2 <- data.table::merge.data.table(map,
                               uberon@elementMetadata,
                               all.x = TRUE,
                               by.x="tissue_ontology_term_id",
                               by.y="id")
  map2[,ancestor:=data.table::fcoalesce(ancestor,tissue_ontology_term_id)]
  map2[,ancestor_name:=data.table::fcoalesce(ancestor_name,as.character(tissue))]
  #### Compute metrics for each celltype ####
  map2[,n_tissues:=data.table::uniqueN(tissue),
      by=c("cell_type_ontology_term_id","cell_type")]
  map2[,n_ancestors:=data.table::uniqueN(ancestor),
       by=c("cell_type_ontology_term_id","cell_type")]




  #### Set cell type order ####
  dend <- KGExplorer::ontology_to(ont = cl,
                                  # terms = results$cell_type_ontology_term_id,
                                  to="dendrogram")
  map2$cell_type_ontology_term_id <- factor(map2$cell_type_ontology_term_id,
                                            levels=labels(dend),
                                            ordered = TRUE)
  #### Create tissue-celltype heatmap #####

  ggplot2::ggplot(map2,
                  ggplot2::aes(x=cell_type, y=ancestor_name,
                               fill=ancestor_name)
                  ) +
    ggplot2::geom_tile(show.legend = FALSE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1, vjust=0.5))


  return(results_cl)
}



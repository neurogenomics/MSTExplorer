#' Get Cell Ontology
#'
#' A thin wrapper around \link[KGExplorer]{get_ontology} to get the
#' specific version of the Cell Ontology used in the original analyses that
#' produced the results stored in \link{load_example_results}.
#'
#' @inheritParams KGExplorer::get_ontology
#' @inheritDotParams KGExplorer::get_ontology
#' @export
#' @examples
#' cl <- get_cl()
get_cl <- function(name = "cl",
                   tag = "v2023-09-21",
                   lvl = 1,
                   remove_rings=TRUE,
                   ...){
  KGExplorer::get_ontology(name = name,
                           tag = tag,
                           lvl = lvl,
                           remove_rings = remove_rings,
                           ...) |>
    KGExplorer::filter_ontology(
      keep_descendants = "cell"
    )
}

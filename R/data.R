#' Example prioritised targets
#'
#' @description
#' Example output from the function \link[MultiEWCE]{prioritise_targets}.
#' Used in examples to reduce run time.
#' @source
#' \code{
#' example_targets <- prioritise_targets()
#' usethis::use_data(example_targets, overwrite = TRUE)
#' }
#' @format ontology_index
#' @usage data("example_targets")
"example_targets"




#' Cell-type maps
#'
#' @description
#' Mappings from cell-type names to cell-type ontology terms in the
#' HumanCellLandscape and DescartesHuman scRNA-seq datasets.
#' @source
#' \code{
#' ## Download standardised Seurat objects from CellXGene.
#' obj_list <- list()
#' obj_list[["DescartesHuman"]] <- readRDS("~/Downloads/db4d63ab-0ca3-4fa5-b0cf-34955227912d.rds")
#' obj_list[["HumanCellLandscape"]] <- readRDS("~/Downloads/1b7484e3-83a0-47fe-847e-54d811a2adae.rds")
#' cols <- c("Main_cluster_name",
#'           "author_cell_type",
#'           "cell_type_ontology_term_id",
#'           "cell_type")
#' celltype_maps <- lapply(obj_list, function(obj){
#'   select_cols <- cols[cols %in% colnames(obj@meta.data)]
#'   meta <- unique(data.table::as.data.table(obj@meta.data[,select_cols]))
#'   data.table::setnames(meta, 1, "author_celltype")
#'   meta[,author_celltype:=EWCE::fix_celltype_names(author_celltype, make_unique = FALSE)]
#'   unique(meta)
#' }) |> data.table::rbindlist(idcol = "ctd")
#' #### Ensure one entry per cell type
#' celltype_maps <- celltype_maps[,.SD[1],by=c("ctd","author_celltype")]
#' usethis::use_data(celltype_maps, overwrite = TRUE)
#' }
#' @format data.table
#' @usage data("celltype_maps")
"celltype_maps"


#' Anatomy maps
#'
#' @description
#' Mappings from cell-type names to anatomy ontology terms in the
#' HumanCellLandscape and DescartesHuman scRNA-seq datasets.
#' @source
#' \code{
#' ## Download standardised Seurat objects from CellXGene.
#' obj_list <- list()
#' obj_list[["DescartesHuman"]] <- readRDS("~/Downloads/db4d63ab-0ca3-4fa5-b0cf-34955227912d.rds")
#' obj_list[["HumanCellLandscape"]] <- readRDS("~/Downloads/1b7484e3-83a0-47fe-847e-54d811a2adae.rds")
#' cols <- c("tissue_ontology_term_id",
#'           "tissue",
#'           # "tissue_original",
#'           "cell_type_ontology_term_id",
#'           "cell_type")
#' anatomy_maps <- lapply(obj_list, function(obj){
#'   select_cols <- cols[cols %in% colnames(obj@meta.data)]
#'   unique(data.table::as.data.table(obj@meta.data[,select_cols]))
#' }) |> data.table::rbindlist(idcol = "ctd", fill=TRUE)
#' usethis::use_data(anatomy_maps, overwrite = TRUE)
#' }
#' @format data.table
#' @usage data("anatomy_maps")
"anatomy_maps"


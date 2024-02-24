#' Terminal cell types
#'
#' A manually-curated list of terminally differentiated cell types
#' along with supporting references from the literature.
#' @param fix_names Standardise celltype names with
#'  \link[EWCE]{fix_celltype_names}.
#' @returns data.table of cell types.
#'
#' @export
#' @importFrom data.table fread
#' @examples
#' celltypes <- terminal_celltypes()
terminal_celltypes <- function(fix_names=TRUE){
  file <- system.file("extdata","terminal_celltypes.csv.gz",
                      package = "MSTExplorer")
  ct <- data.table::fread(file)
  if(isTRUE(fix_names)){
    ct$CellType <- EWCE::fix_celltype_names(celltypes = ct$CellType,
                                            make_unique = FALSE)
  }
  return(ct)
}

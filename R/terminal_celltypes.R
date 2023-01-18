#' Terminal cell types
#'
#' A manually-curated list of terminally differentiated cell types
#' along with supporting references from the literature.
#' @returns data.table of cell types.
#'
#' @export
#' @importFrom data.table fread
#' @examples
#' celltypes <- terminal_celltypes()
terminal_celltypes <- function(){
  file <- system.file("extdata","terminal_celltypes.csv.gz",
                      package = "MultiEWCE")
  data.table::fread(file)
}

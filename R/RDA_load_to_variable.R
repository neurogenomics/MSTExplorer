#' Load RDA file and assign to specific variable
#' Can be useful, but sort of pointless now as we have decided to only use .rds
#'
#' @param file file path to rda \<string\>
#' @examples
#'
#' \dontrun{
#' files = "results.rds"
#' RDA_assign_load(file)}
#'
#' @returns the data contained in the rda file
#' @export
RDA_assign_load <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  return(tmp[[ls(tmp)[1]]])
}

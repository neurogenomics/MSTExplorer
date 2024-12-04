#' Extract help
#'
#' Extract help documentation from a function.
#' @source \href{https://stackoverflow.com/a/9195691}{Original source.}
#' @param pkg Package name.
#' @param fn Function name,
#' @returns text. html, latex formats
#'
#' @keywords internal
#' @importFrom utils getFromNamespace capture.output
extract_help <- function(pkg,
                         fn = NULL,
                         to = c("txt", "html", "latex", "ex")){
  requireNamespace("tools")

  to <- match.arg(to)
  rdbfile <- file.path(find.package(pkg, lib.loc = .libPaths()), "help", pkg)
  fetchRdDB <- utils::getFromNamespace("fetchRdDB","tools")
  rdb <- fetchRdDB(rdbfile, key = fn)
  convertor <- switch(to,
                      txt   = tools::Rd2txt,
                      html  = tools::Rd2HTML,
                      latex = tools::Rd2latex,
                      ex    = tools::Rd2ex
  )
  f <- function(x) utils::capture.output(convertor(x))
  if(is.null(fn)) lapply(rdb, f) else f(rdb)
}

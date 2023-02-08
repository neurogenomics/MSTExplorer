extract_help <- function(pkg,
                         fn = NULL,
                         to = c("txt", "html", "latex", "ex")){
  # from: https://stackoverflow.com/a/9195691
  to <- match.arg(to)
  rdbfile <- file.path(find.package(pkg), "help", pkg)
  rdb <- tools:::fetchRdDB(rdbfile, key = fn)
  convertor <- switch(to,
                      txt   = tools::Rd2txt,
                      html  = tools::Rd2HTML,
                      latex = tools::Rd2latex,
                      ex    = tools::Rd2ex
  )
  f <- function(x) capture.output(convertor(x))
  if(is.null(fn)) lapply(rdb, f) else f(rdb)
}

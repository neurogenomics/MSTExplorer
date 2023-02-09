extract_filters <- function(pkg = "MultiEWCE",
                            fn = "prioritise_targets"){

  requireNamespace("rvest")
  level <- NULL;

  out <- extract_help(pkg = pkg,
                      fn = fn,
                      to = "html")
  filters <- (
    rvest::read_html(paste(out,collapse = "")) |>
      rvest::html_elements("ol")
  )[[1]] |>
    rvest::html_children() |>
    rvest::html_text() |>
    stringr::str_split(":",n = 3, simplify = TRUE) |>
    data.table::as.data.table() |>
    `colnames<-`(c("level","step","description"))
  cols <- names(filters)
  filters[ , (cols) := lapply(.SD, trimws), .SDcols = cols]
  filters[,level:=gsub("[ ]+|[-]|level|levels","",level)]
  #### Print for plot legend ####
  # cat(paste(filters[,legend:=paste0(step,": ",description)]$legend, collapse = "\n"))
  return(filters)
}

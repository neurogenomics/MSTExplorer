extract_filters <- function(pkg = "MSTExplorer",
                            fn = "prioritise_targets"){

  requireNamespace("rvest")
  level <- NULL;

  out <- extract_help(pkg = pkg,
                      fn = fn,
                      to = "html")

  txt <- (rvest::read_html(paste(out,collapse = ""))|>
    rvest::html_elements("dl"))[[1]] |>
    rvest::html_children()|>
    rvest::html_text()

  filters <- data.table::data.table(
    stringr::str_split(txt[seq(from = 1, to = length(txt)-1, by = 2)] ,"; ",
                       n = 2,
                       simplify = TRUE)[,c(1,2)],
    txt[seq(from = 2, to = length(txt), by = 2)]
  )|>
    `colnames<-`(c("level","step","description"))

  cols <- names(filters)
  filters[ , (cols) := lapply(.SD, trimws), .SDcols = cols]
  filters[ , (cols) := lapply(.SD, trimws, whitespace=";"), .SDcols = cols]
  filters[,level:=gsub("[-]|level|levels","",level)]
  #### Print for plot legend ####
  # cat(paste(filters[,legend:=paste0(step,": ",description)]$legend, collapse = "\n"))
  return(filters)
}

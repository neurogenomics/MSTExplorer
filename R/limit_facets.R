limit_facets <- function(dat,
                         facet_var,
                         facet_subset=NULL,
                         max_facets = NULL){
  if(!is.null(facet_subset)){
    dat <- dat[get(facet_var) %in% facet_subset,]
  }
  if(!is.null(max_facets)){
    ids <- utils::head(unique(dat[[facet_var]]),max_facets)
    dat <- dat[get(facet_var) %in% ids,]
  }
  return(dat)
}

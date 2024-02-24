construct_labeller <- function(dat,
                               facets,
                               facets_n,
                               suffix="phenotypes"){
  facet_label <- NULL;
  if(is.null(facets_n) ||
     !facets_n %in% names(dat)) {
    "label_value"
  } else {
    dat[,facet_label:=paste0(
      get(facets)," (n=",get(facets_n),
      if(is.null(suffix))"" else paste0(" ",suffix),")")
    ]
    ggplot2::labeller(
      .cols = Reduce(
        stats::setNames,
        as.list(unique(dat[,c("facet_label",facets),with=FALSE]))
        )
    )
  }
}

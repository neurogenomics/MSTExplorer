make_rannot <- function(annot,
                        row_side_vars,
                        palettes = list(pals::kovesi.cyclic_mrybm_35_75_c68_s25,
                                        pals::kovesi.linear_bgy_10_95_c74,
                                        pals::kovesi.linear_bmw_5_95_c86)
                        ){
  annot <- annot[,row_side_vars,with=FALSE]
  cat_cols <- names(annot)#[!unlist(lapply(annot, is.numeric))]
  col_annot <- lapply(stats::setNames(cat_cols,
                                      cat_cols), function(x){
                i <- which(names(annot)==x)
                vals <- as.character(sort(unique(annot[,i,with=FALSE][[1]])))
                dict <- stats::setNames(palettes[[i]](n = length(vals)),
                                        vals)
                dict[vals]
              })
  ra <- ComplexHeatmap::HeatmapAnnotation(
    df = annot,
    which = "row",
    col =  col_annot)
  return(ra)
}

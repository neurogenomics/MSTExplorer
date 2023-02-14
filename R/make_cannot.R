make_cannot <- function(annot,
                        col_side_vars,
                        palettes = list(pals::kovesi.cyclic_mrybm_35_75_c68_s25,
                                        pals::kovesi.linear_bgy_10_95_c74,
                                        pals::kovesi.linear_bmw_5_95_c86)
                        ){
  ha_list <- lapply(stats::setNames(col_side_vars,
                                    gsub("celltype","cell type",
                                         gsub("_"," ",col_side_vars)
                                    )
  ), function(x){
    i <- which(col_side_vars==x)
    vals <- annot[[x]]
    dict <- stats::setNames(palettes[[i]](n = length(vals)),
                            vals)
    ComplexHeatmap::anno_barplot(which = "column",
                                 x = vals,
                                 gp = grid::gpar(fill = dict[vals]),
                                 add_numbers = FALSE)
  } )
  ha <- do.call(ComplexHeatmap::columnAnnotation, ha_list )
  # ComplexHeatmap::draw(ha)
  return(ha)
}

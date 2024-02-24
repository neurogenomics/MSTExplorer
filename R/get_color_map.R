get_color_map <- function(dat,
                          columns = "top_ancestor_name",
                          ddata,
                          celltype_col,
                          celltype_col_order=ddata$labels$label,
                          preferred_palettes="tol"){
  ancestor_color <- top_ancestor_name <- NULL;
  color_map <- KGExplorer::map_colors(dat,
                                      columns = columns,
                                      preferred_palettes = preferred_palettes,
                                      as="dict")[[1]]
  ####sort dat rows by levels in ddata$labels$id
  dat2 <- unique(dat[,c(celltype_col,columns), with=FALSE])
  if(!is.null(celltype_col_order)){
    dat2 <- dat2[order(match(dat2[[celltype_col]],celltype_col_order)),]
  }
  dat2[,ancestor_color:=color_map[top_ancestor_name]]
  color_vector <- dat2$ancestor_color
  return(list(
    color_map=color_map,
    color_vector=color_vector
  ))
}

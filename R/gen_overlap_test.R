gen_overlap_test <- function(ct_genes,
                             dgenes,
                             list_name=NULL,
                             long_format=FALSE,
                             bg,
                             verbose=TRUE){
  if(!is.null(list_name)){
    messager("Testing overlap: ",unique(list_name), parallel = TRUE,
             v = verbose)
  }
    lapply(ct_genes, function(cgenes){
      r <- GeneOverlap::newGeneOverlap(
        listA = cgenes,
        listB = dgenes,
        genome.size = length(bg)) |>
        GeneOverlap::testGeneOverlap()
      data.table::data.table(
        intersection=if(long_format) r@intersection else list(r@intersection),
        intersection_size=length(r@intersection),
        # union=if(long_format) r@union else list(r@union),
        union_size=length(r@union),
        pval=r@pval,
        odds.ratio=r@odds.ratio,
        jaccard=r@Jaccard)
    }) |>
      data.table::rbindlist(use.names = TRUE, idcol = "CellType")
}

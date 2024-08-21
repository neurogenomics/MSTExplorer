prepare_congenital_annotations <- function(results,
                                           fetal_keywords,
                                           celltype_col="author_celltype",
                                           gpt_annot = HPOExplorer::gpt_annot_codify()){
  fetal_celltype <- fetal_nonfetal_pdiff <- stage <-
    has_adult_and_fetal <- NULL;
  #### Prepare data ####
  results <- HPOExplorer::add_gpt_annotations(results,
                                              annot = gpt_annot$annot)
  results <- map_celltype(results)
  if(!is.null(fetal_keywords)){
    results[,fetal_celltype:=grepl(paste(fetal_keywords,collapse="|"),
                                   get(celltype_col),
                                   ignore.case = TRUE)]
  } else {
    results[,fetal_celltype:=stage!="Adult"]
  }
  results[,has_adult_and_fetal:=(
    TRUE %in% fetal_celltype & FALSE %in% fetal_celltype
  ), by=c("hpo_id","cl_name")]
  #### Compute difference between fetal/nonfetal cell types ####
  results[has_adult_and_fetal==TRUE,
          fetal_nonfetal_pdiff:=(
            mean(p[fetal_celltype==FALSE], na.rm=TRUE) -
              mean(p[fetal_celltype==TRUE], na.rm=TRUE)
          ), by=c("hpo_id","cl_id")]
  return(results)
}

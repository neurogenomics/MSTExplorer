#' Get target branches
#'
#' A named list of mappings between HPO ancestral terms (list names),
#' and Cell Ontology ancestral terms (list values).
#' @keywords internal
get_target_branches <- function(sort=TRUE){
  target_branches <- list(
    "Abnormality of the nervous system"=
      c("neural cell"),
    "Abnormality of the cardiovascular system"=
      c("cardiocyte"),
    "Abnormality of the immune system"=
      c("leukocyte"),
    "Abnormality of the musculoskeletal system"=
      c("cell of skeletal muscle",
        # "bone cell",
        "chondrocyte"
        ),
    "Abnormality of the respiratory system"=
      c(
        # "ciliated epithelial cell",
        "respiratory epithelial cell",
        "epithelial cell of lung"
        ),
    "Abnormality of the endocrine system"=
      c("endocrine cell"),
    "Abnormality of the eye"=c(
      "photoreceptor cell",
      "retinal cell"
    )
  )
  if(isTRUE(sort)){
    target_branches <- target_branches[sort(names(target_branches))]
  }
  return(target_branches)
}

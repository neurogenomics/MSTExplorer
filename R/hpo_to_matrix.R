#' HPO to matrix
#'
#' Convert gene-phenotype associations from the Human Phenotype Ontology (HPO)
#' into a gene x phenotype matrix.
#' @param terms A subset of HPO IDs to include.
#' Set to \code{NULL} (default) to include all terms.
#' @param run_cor Return a matrix of pairwise correlations.
#' @inheritParams HPOExplorer::make_phenos_dataframe
#' @inheritParams data.table::dcast.data.table
#' @inheritParams stats::cor
#' @returns A gene x phenotype matrix,
#' or a phenotype x phenotype matrix if \code{run_cor=TRUE}.
#'
#' @keywords internal
#' @importFrom HPOExplorer load_phenotype_to_genes
#' @importFrom data.table dcast.data.table copy setnafill :=
#' @importFrom stats terms as.formula cor
hpo_to_matrix <- function(terms = NULL,
                          phenotype_to_genes =
                            HPOExplorer::load_phenotype_to_genes(),
                          formula = "Gene ~ Phenotype",
                          fun.aggregate = mean,
                          fill = 0,
                          run_cor = FALSE,
                          as_matrix = TRUE,
                          sparse = FALSE,
                          method = "pearson",
                          verbose = TRUE){

  ID <- dummy <- NULL;

  messager("Constructing HPO gene x phenotype matrix.",v=verbose)
  if(!is.null(terms)){
    phenotype_to_genes <- phenotype_to_genes[ID %in% unique(terms),]
  }
  #### Cast into gene x phenotype matrix ####
  X_dt <- phenotype_to_genes[,dummy:=1] |>
    data.table::dcast.data.table(formula = formula,
                                 value.var = "dummy",
                                 fun.aggregate = fun.aggregate,
                                 fill = fill,
                                 na.rm = TRUE)
  #### Make matrix nownames ####
  meta_vars <- all.vars(stats::terms(stats::as.formula(formula))[-1])
  rn <- data.table::copy(X_dt)[, rn:=do.call(paste0,.SD),
                               .SDcols=meta_vars]$rn
  #### Fill NAs ####
  if(!is.null(fill)){
    data.table::setnafill(X_dt, fill = fill,
                          cols = names(X_dt[,-meta_vars, with=FALSE]))
  }
  #### Format and return ####
  if(isTRUE(as_matrix)){
    X <- as.matrix(X_dt[,-meta_vars,with=FALSE])|> `rownames<-`(rn)
    if(isTRUE(run_cor)){
      messager("Computing all parwise correlations.",v=verbose)
      X_cor <- stats::cor(X, method = method)
      # stats::heatmap(X_cor)
      return(X_cor)
    } else {
      return(X)
    }
  } else{
    return(X_dt)
  }
}

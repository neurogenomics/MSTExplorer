run_similarity <- function(X_list=NULL,
                           method="spearman",
                           nThreads=0){

  Xcor_list <- list()
  if(is.null(X_list)){
    messager("Creating default correlation matrices.")
    X_list <- list()
    #### Ontology-based similarity ####
    X_list[["ontology"]] <- HPOExplorer::get_hpo()
    Xcor_list[["ontology"]] <- KGExplorer::ontology_to(ont=X_list[["ontology"]],
                                                       to="similarity")
    #### Gene-based similarity ####
    X_list[["genes"]] <- HPOExplorer::hpo_to_matrix(
      formula = "gene_symbol ~ hpo_id"
      )
    Xcor_list[["genes"]] <- WGCNA::cor(X_list[["genes"]],
                                       nThreads=nThreads,
                                       method=method)
    #### Celltype-based similarity ####
    results <- load_example_results()
    # results[,combined_score:=(1-p)*fold_change]
    X_list[["celltypes"]] <- data.table::dcast.data.table(
      results,
      formula =  paste0(ctd,".",CellType) ~ hpo_id,
      value.var = "p",
      fill = 1,
      fun.aggregate = mean,
      na.rm=TRUE
      ) |>
      as.data.frame() |>
      tibble::column_to_rownames("ctd") |>
      as.matrix()|>
      methods::as("sparseMatrix")
    Xcor_list[["celltypes"]] <- WGCNA::cor(X_list[["celltypes"]],
                                           nThreads=nThreads,
                                           method=method)
  }
  #### Find intersecting rownames ####
  ids <- Reduce(intersect, lapply(Xcor_list, rownames))
  ### Test how similar each pair of matrices are to each other ###
  cor_tests <- combn(x=names(Xcor_list), m=2, simplify=FALSE,
                        function(x){
    messager("Testing:", x[1], "vs.", x[2])
      stats::cor.test(Xcor_list[[x[1]]][ids, ids],
                      Xcor_list[[x[2]]][ids, ids],
             method="kendall")
  })|> `names<-`( combn(x=names(Xcor_list), m=2, simplify=FALSE,
                        function(x) paste(x, collapse = ".")) )
  cor_dt <- lapply(cor_tests, broom::tidy)|>
    data.table::rbindlist(idcol = "test")



  return(
    list(X_list = X_list,
         Xcor_list = Xcor_list,
         cor_dt = cor_dt)
  )
}

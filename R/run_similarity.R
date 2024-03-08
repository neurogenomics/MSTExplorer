run_similarity <- function(X_list=NULL){
  if(is.null(X_list)){
    messager("Creating default correlation matrices.")
    X_list <- list()
    #### Ontology-based similarity: ground-truth phenotype similarity ####
    hpo <- HPOExplorer::get_hpo()
    X_list[["ontology"]] <- KGExplorer::ontology_to(ont=hpo,
                                                    to="similarity")
    #### Gene-based similarity ####
    X_list[["genes"]] <- HPOExplorer::hpo_to_matrix(
      formula = "gene_symbol ~ hpo_id"
      )
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
  }
  #### Run Correlations on each matrix ####
  Xcor_list <- lapply(X_list, function(X) {
    WGCNA::cor(X)
  })
  #### Run decomposition before computing correlations ####
  reduc_list <- lapply(X_list, function(X) {
    # fastICA::fastICA(Matrix::t(X),
    #                  method="C",
    #                  n.comp=100)
    # a$X %*% a$K
    phenomix::run_pca(mat = X,
                      transpose = TRUE,
                      ncomp = 100)
  })
  # tmp <- Seurat::RunPCA( Xcor_list[["ontology"]])
  # pca_dt <- data.frame(variance_explained=(pca$sdev^2/sum(pca$sdev^2))[1:100],
  #                      PC=1:100) |>
  #   dplyr::mutate(type="genes")
  # ggplot2::ggplot(pca_dt, ggplot2::aes(x=PC, y=variance_explained, color=type)) +
  #   ggplot2::geom_line() +
  #   ggplot2::theme_minimal() +
  #   ggplot2::labs(title="Variance explained by principal components",
  #                 x="Principal component",
  #                 y="Variance explained") +
  #   ggplot2::scale_color_manual(values=c("blue","red","green"))


  #### Find intersecting rownames ####
  ids <- Reduce(intersect, lapply(Xcor_list, rownames))
  #### Test how similar each pair of correlation matrices are to each other ####
  cor_tests <- combn(x=names(Xcor_list), m=2, simplify=FALSE,
                        function(x){
    messager("Testing:", x[1], "vs.", x[2])
      stats::cor.test(Xcor_list[[x[1]]][ids, ids],
                      Xcor_list[[x[2]]][ids, ids],
             method="pearson")
  })|> `names<-`( combn(x=names(Xcor_list), m=2, simplify=FALSE,
                        function(x) paste(x, collapse = ".")) )
  cor_dt <- lapply(cor_tests, broom::tidy)|>
    data.table::rbindlist(idcol = "test")


  #### Return ####
  return(
    list(X_list = X_list,
         Xcor_list = Xcor_list,
         cor_dt = cor_dt)
  )
}

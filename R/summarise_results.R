#' Summarise results
#'
#' Summarise results from \pkg{MSTExplorer}.
#' @export
summarise_results <- function(results,
                              group_var="ctd",
                              digits=3,
                              save_path=tempfile("summarise_results.csv")){

  CellType <- hpo_id <- celltypes_per_phenotype <- phenotypes_per_celltype <- NULL;
  p2g <- HPOExplorer::load_phenotype_to_genes()
  p2g[,phenotype:=hpo_id]
  total_phenotypes <- data.table::uniqueN(p2g$hpo_id)
  total_diseases <- data.table::uniqueN(p2g$disease_id)
  diseases_covered <- data.table::uniqueN(p2g[hpo_id %in% results[q<0.05]$hpo_id]$disease_id)
  #### table 1 ####
  t1 <- results[,list(
    ## tests
    `tests significant`=sum(q<0.05),
    `tests`=.N,
    `tests significant (%)`=100*sum(q<0.05)/.N,
    ## celltypes
    `cell types significant`=data.table::uniqueN(CellType[q<0.05]),
    `cell types`=data.table::uniqueN(CellType),
    `cell types significant (%)`=
      100*data.table::uniqueN(CellType[q<0.05])/data.table::uniqueN(CellType),
    ## phenotypes
    `phenotypes significant`=data.table::uniqueN(hpo_id[q<0.05]),
    `phenotypes tested`=data.table::uniqueN(hpo_id),
    `phenotypes`=total_phenotypes,
    `phenotypes significant (%)`=
      100*data.table::uniqueN(hpo_id[q<0.05])/total_phenotypes
  ),
  by=group_var]|>unique()
  #### table 2 ####
  t2 <- results[,

          list(
            ## diseases
            `diseases significant`= p2g[phenotype %in% hpo_id[q<0.05]]$disease_id |> data.table::uniqueN(),
            `diseases`=p2g$disease_id|>data.table::uniqueN(),
            `diseases significant (%)`=100*p2g[phenotype %in% hpo_id[q<0.05]]$disease_id |> data.table::uniqueN()/total_diseases
          ), by=group_var
  ]|> unique()
  #### table 3 ####
  t3 <- results[,
          list(
            `celltypes_per_phenotype`=data.table::uniqueN(CellType[q<0.05])
          ),
          by=c(group_var,"hpo_id")
          ][,
            list(
              `cell types per phenotype (mean)`=mean(celltypes_per_phenotype),
              `cell types per phenotype (median)`=median(celltypes_per_phenotype),
              `cell types per phenotype (min)`=min(celltypes_per_phenotype),
              `cell types per phenotype (max)`=max(celltypes_per_phenotype)
            ),
            by=group_var
          ]
  #### table 3 ####
  t4 <- results[,
                list(
                  `phenotypes_per_celltype`=data.table::uniqueN(hpo_id[q<0.05])
                ),
                by=c(group_var,"CellType")
  ][,
    list(
      `phenotypes per cell type (mean)`=mean(phenotypes_per_celltype),
      `phenotypes per cell type (median)`=median(phenotypes_per_celltype),
      `phenotypes per cell type (min)`=min(phenotypes_per_celltype),
      `phenotypes per cell type (max)`=max(phenotypes_per_celltype)
    ),
    by=group_var
  ]

  #### merge tables ####
  tmerged <- merge(t1,t2, by=group_var)|>
    merge(t3, by=group_var) |>
    merge(t4, by=group_var) |>
    format(big.mark=",", digits=digits)
  tmerged_transposed <- t(tmerged)
  #### Save ####
  if(!is.null(save_path)){
    utils::write.csv(tmerged_transposed,"~/Downloads/tmerged.csv")
  }
  #### Return ####
  return(list(
    tmerged=data.table::data.table(tmerged),
    tmerged_transposed=tmerged_transposed
  ))
}
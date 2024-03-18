#' Validate associations with the Monarch Knowledge Graph
#'
#' This function validates the associations between phenotypes and cell types
#' with the Monarch Knowledge Graph (MKG). It computes the number of phenotypes
#' and cell types that are present in the MKG and the number of associations
#' that are present in the MKG.
#' @inheritParams prioritise_targets
#' @param kg A data.table with phenotype-cell type relationships
#' ("from" and "to" columns, respectively) gathered from
#'  the Monarch Knowledge Graph.
#' @export
#' @examples
#' results <- load_example_results()
#' kg <- validate_associations_mkg(result=results)
validate_associations_mkg <- function(results=load_example_results(),
                                      kg=get_data("monarch_kg_cells.csv"),
                                      q_threshold=0.05){
  from <- hpo_id <- NULL;
  # kg <- data.table::fread(here::here("data/monarch_kg_cells.csv"))
  add_logfc(results)
  kg_og <- data.table::copy(kg)
  kg <- kg[grepl("HP:",from)][from %in% unique(results$hpo_id)]
  message(paste(
    "Remaining:",length(unique(kg$from)),
    "phenotypes across",length(unique(kg$to)),"celltypes."
  ))
  #### Merge results with Knowledge Graph ####
  kg_res <- data.table::merge.data.table(
    kg[,c("from","to","label.to")],
    results[q<q_threshold],
    by.x="from",
    by.y="hpo_id") |>
    data.table::setorderv(c("from","logFC"), c(1,-1))
  #### Compute average number of celltypes / phenotype ####
  # kg[,list(N=data.table::uniqueN(to)), by="from"]$N|>summary()
  # kg_res[,list(N=data.table::uniqueN(CellType)), by="from"]$N|>summary()

  missing_phenos <- setdiff(kg$from, kg_res$from)
  proportion_phenotypes_captured <- (1-(length(missing_phenos) /
                                        length(unique(kg$from))))
  message(format(proportion_phenotypes_captured*100,digits = 4),
          "% phenotypes recovered.")
  #### Get missing phenotypes and associations ####
  kg_missing <- kg[from%in% missing_phenos]
  res_missing <- results[hpo_id%in% missing_phenos]
  #### Compute increase in knowledge ####
  # n_before <- kg$from|>unique()|>length()
  # n_after <- res[q<0.05]$hpo_id|>unique()|>length()
  n_before <- nrow(kg)
  n_after <- nrow(results[q<q_threshold])

  return(
    list(
      kg=kg_og,
      kg_filt=kg,
      kg_res=kg_res,
      proportion_phenotypes_captured=proportion_phenotypes_captured,
      missing_phenos=missing_phenos,
      kg_missing=kg_missing,
      res_missing=res_missing,
      n_before=n_before,
      n_after=n_after,
      increased_knowledge=n_after/n_before
    )
  )
}

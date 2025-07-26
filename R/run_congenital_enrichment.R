#' Run congenital enrichment analyses
#'
#' First, this function computes the mean difference in association p-values
#' between a given phenotype and the congenital and non-congenital versions
#' of a cell type (the \code{fetal_nonfetal_pdiff} column).
#' It then sorts the results from the largest to the smallest difference.
#' Here, a positive difference indicates that the phenotype is more associated
#' with the foetal/embryonic version of the cell type, while a negative difference
#' indicates that the phenotype is more associated with the adult version.
#' Next, it runs enrichment analyses for the top and bottom N phenotypes
#' using the \link[simona]{dag_enrich_on_offsprings} function.
#' Finally, it run enrichment analyses for the top and bottom N cell types.
#' @param top_n_phenotypes Number of top and bottom phenotypes to return.
#' @param top_n_hpo Number of top and bottom phenotypes
#'  to run enrichment analyses on.
#' @param top_n_cl Number of top and bottom cell types
#'  to run enrichment analyses on.
#' @param prune Prune redundant ancestral terms from the enrichment results
#'  using \link[KGExplorer]{prune_ancestors}.
#' @inheritParams prioritise_targets
#' @inheritParams plot_congenital_annotations
#' @inheritParams plot_bar_dendro
#' @inheritDotParams simona::dag_enrich_on_offsprings
#' @returns Named list of enrichment results.
#'
#' @export
#' @examples
#' results <- load_example_results()[ctd=="HumanCellLandscape"]
#' out <- run_congenital_enrichment(results=results)
run_congenital_enrichment <- function(results,
                                      hpo=HPOExplorer::get_hpo(),
                                      cl=get_cl(),
                                      gpt_annot = HPOExplorer::gpt_annot_codify(),
                                      fetal_keywords=c("fetal",
                                                       "fetus",
                                                       "primordial",
                                                       "hESC",
                                                       "embryonic"),
                                      celltype_col="author_celltype",
                                      top_n_phenotypes=10,
                                      top_n_hpo=50,
                                      top_n_cl=top_n_hpo,
                                      q_threshold=.05,
                                      prune=TRUE,
                                      ...){

  p_adjust <- has_adult_and_fetal <- fetal_nonfetal_pdiff <- NULL;

  results <- prepare_congenital_annotations(results = results,
                                            fetal_keywords = fetal_keywords,
                                            celltype_col = celltype_col,
                                            gpt_annot = gpt_annot)
  fetal_pdiff_top <-
    results[fetal_nonfetal_pdiff >= min(utils::tail(sort(unique(results$fetal_nonfetal_pdiff)), top_n_phenotypes))]|>
    data.table::setorderv("fetal_nonfetal_pdiff", -1)
  fetal_pdiff_bottom <-
    results[fetal_nonfetal_pdiff <= max(utils::head(sort(unique(results$fetal_nonfetal_pdiff)), top_n_phenotypes))]|>
    data.table::setorderv("fetal_nonfetal_pdiff", 1)

  gdat <- unique(results[has_adult_and_fetal==TRUE,
                         c("hpo_id","cl_id","cl_name",
                           "fetal_nonfetal_pdiff",
                           "congenital_onset")
  ][!is.na(fetal_nonfetal_pdiff)])|>
    unique() |>
    data.table::setorderv("fetal_nonfetal_pdiff", -1)
  gdat_by_pheno <- results[,list(fetal_nonfetal_pdiff=mean(fetal_nonfetal_pdiff, na.rm=TRUE)),
                           by=c("hpo_id","congenital_onset")]|>
    data.table::setorderv("fetal_nonfetal_pdiff", -1)
  gdat_by_celltype <- results[,list(fetal_nonfetal_pdiff=mean(fetal_nonfetal_pdiff, na.rm=TRUE)),
                              by=c("cl_name","cl_id")]|>
    data.table::setorderv("fetal_nonfetal_pdiff", -1)
  #### HPO enrichments ####
  hpo_enrich_top <- (
    simona::dag_enrich_on_offsprings(
      dag=hpo,
      terms = utils::head((gdat$hpo_id),top_n_hpo),
      ...)|>
        data.table::data.table(key = "p_value")
    )[p_adjust<q_threshold]
  hpo_enrich_bottom <- (
    simona::dag_enrich_on_offsprings(
      dag=hpo,
      terms = utils::tail((gdat$hpo_id),top_n_hpo),
      ...)|>
        data.table::data.table(key = "p_value")
    )[p_adjust<q_threshold]
  #### CL enrichments ####
  # cl2 <- KGExplorer::filter_ontology(cl, terms = as.character(unique(gdat$cl_id)))
  cl_enrich_top <- (
    simona::dag_enrich_on_offsprings(
      dag=cl,
      terms = as.character(utils::head((gdat$cl_id),top_n_cl)),
      ...)|>
        data.table::data.table(key = "p_value")
    )[p_adjust<q_threshold]
  cl_enrich_bottom <- (
    simona::dag_enrich_on_offsprings(
      dag=cl,
      terms = as.character(utils::tail((gdat_by_celltype$cl_id),top_n_cl)),
      ...)|>
        data.table::data.table(key = "p_value")
    )[p_adjust<q_threshold]

  if(isTRUE(prune)){
    hpo_enrich_top = KGExplorer::prune_ancestors(hpo_enrich_top,
                                                 ont=hpo,
                                                 id_col="term")
    hpo_enrich_bottom = KGExplorer::prune_ancestors(hpo_enrich_bottom,
                                                 ont=hpo,
                                                 id_col="term")
    cl_enrich_top = KGExplorer::prune_ancestors(cl_enrich_top,
                                                 ont=cl,
                                                 id_col="term")
    cl_enrich_bottom = KGExplorer::prune_ancestors(cl_enrich_bottom,
                                                    ont=cl,
                                                    id_col="term")
  }

  #### GLM test ####
  #### Test for difference between congenital frequency of phenotype ####
  ##  linear (.L), the second is quadratic (.Q), the third is cubic (.C), the fourth is quartic (Year^4),
  # gres <- lme4::lmer(data=gdat, fetal_nonfetal_pdiff ~ congenital_onset | cl_id)
  # gres <- lm(data=gdat, fetal_nonfetal_pdiff ~ congenital_onset)
  # # broom::tidy(gres)
  # summary(gres)
  #
  # ggplot2::ggplot(gdat, ggplot2::aes(x=congenital_onset, y=fetal_nonfetal_pdiff)) +
  #   ggplot2::geom_boxplot()

  return(
    list(
      fetal_pdiff_top = fetal_pdiff_top,
      fetal_pdiff_bottom = fetal_pdiff_bottom,
      hpo_enrich_top = hpo_enrich_top,
      hpo_enrich_bottom = hpo_enrich_bottom,
      cl_enrich_top = cl_enrich_top,
      cl_enrich_bottom = cl_enrich_bottom
    )
  )
}


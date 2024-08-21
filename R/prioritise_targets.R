#' Prioritise target genes
#'
#' Prioritise target genes based on a procedure:\cr
#' \enumerate{
#' \item{Disease-level: \code{keep_deaths}: }{
#' Keep only diseases with a certain age of death.}
#' \item{Disease-level: \code{severity_threshold_max}: }{
#' \preformatted{
#' Keep only diseases annotated as a certain degree of severity or greater
#'  (filters on maximum severity per disease).}
#' }
#' \item{Phenotype-level: \code{prune_ancestors}: }{
#' \preformatted{
#' Remove redundant ancestral phenotypes when at least one of their
#'  descendants already exist.}
#' }
#' \item{Phenotype-level: \code{keep_descendants}: }{
#' \preformatted{
#' Remove phenotypes belonging to a certain branch of the HPO,
#'  as defined by an ancestor term.}
#' }
#' \item{Phenotype-level: \code{keep_ont_levels}: }{
#' Keep only phenotypes at certain absolute ontology levels within the HPO.}
#' \item{Phenotype-level: \code{pheno_ndiseases_threshold}: }{
#' The maximum number of diseases each phenotype can be associated with.}
#' \item{Phenotype-level: \code{keep_tiers}: }{
#' Keep only phenotypes with high severity Tiers.}
#' \item{Phenotype-level: \code{severity_threshold}: }{
#' Keep only phenotypes with mean Severity equal to or below the threshold.}
#' \item{Phenotype-level: \code{gpt_filters}: }{
#' \preformatted{
#' Keep only phenotypes with certain GPT annotations in specific
#'  severity metrics.}
#' }
#' \item{Phenotype-level: \code{severity_score_gpt_threshold}: }{
#' Keep only phenotypes with a minimum GPT severity score.}
#' \item{Phenotype-level: \code{info_content_threshold}: }{
#' \preformatted{
#' Keep only phenotypes with a minimum information criterion score
#'  (computed from the HPO).}
#' }
#' \item{Symptom-level: \code{pheno_frequency_threshold}: }{
#' \preformatted{
#' Keep only phenotypes with mean frequency equal to or above the threshold
#'  (i.e. how frequently a phenotype is associated with any diseases in
#'  which it occurs).}
#'  }
#' \item{Symptom-level: \code{keep_onsets}: }{
#' Keep only symptoms with a certain age of onset.}
#' \item{Symptom-level: \code{symptom_p_threshold}: }{
#' Uncorrected p-value threshold to filter cell type-symptom associations by.}
#' \item{Symptom-level: \code{symptom_intersection_threshold}: }{
#' \preformatted{
#' Minimum proportion of genes overlapping between a symptom gene list
#'  (phenotype-associated genes in the context of a particular disease)
#'  and the phenotype-cell type association driver genes.}
#'  }
#' \item{Cell type-level: \code{q_threshold}: }{
#' \preformatted{
#' Keep only cell type-phenotype association results at q<=0.05.}
#' }
#' \item{Cell type-level: \code{effect_threshold}: }{
#' Keep only cell type-phenotype association results at effect size>=1.}
#' \item{Cell type-level: \code{keep_celltypes}: }{
#' Keep only terminally differentiated cell types.}
#' \item{Gene-level: \code{keep_chr}: }{
#' Remove genes on non-standard chromosomes.}
#' \item{Gene-level: \code{evidence_score_threshold}: }{
#' \preformatted{
#' Remove genes that are below an aggregate phenotype-gene
#'  evidence score threshold.}
#' }
#' \item{Gene-level: \code{gene_size}: }{
#' Keep only genes <4.3kb in length.}
#' \item{Gene-level: \code{add_driver_genes}: }{
#' \preformatted{
#' Keep only genes that are driving the association with a given phenotype
#'  (inferred by the intersection of phenotype-associated genes and gene with
#'  high-specificity quantiles in the target cell type).}
#' }
#' \item{Gene-level: \code{keep_biotypes}: }{
#' Keep only genes belonging to certain biotypes.}
#' \item{Gene-level: \code{gene_frequency_threshold}: }{
#' \preformatted{
#' Keep only genes at or above a certain mean frequency threshold
#'  (i.e. how frequently a gene is associated with a given phenotype
#'  when observed within a disease).}
#' }
#' \item{Gene-level: \code{keep_specificity_quantiles}: }{
#' \preformatted{
#' Keep only genes in top specificity quantiles
#'  from the cell type dataset (CTD).}
#' }
#' \item{Gene-level: \code{keep_mean_exp_quantiles}: }{
#' \preformatted{
#' Keep only genes in top mean expression quantiles
#'  from the cell type dataset (CTD).}
#' }
#' \item{Gene-level: \code{symptom_gene_overlap}: }{
#' \preformatted{
#' Ensure that genes nominated at the phenotype-level also
#'  appear in the genes overlapping at the cell type-specific symptom-level.}
#' }
#' \item{All levels: \code{sort_cols}: }{
#' \preformatted{
#' Sort candidate targets by one or more columns
#'  (e.g. "severity_score_gpt", "q").
#' }
#' }
#' \item{All levels: \code{top_n}: }{
#' \preformatted{
#' Only return the top N targets per variable group
#'  (specified with the "group_vars" argument).
#'  For example, setting "group_vars" to "hpo_id" and "top_n" to 1 would
#'  only return one target (row) per phenotype ID after sorting.}
#' }
#' }
#'
#' Term key:\cr
#' \itemize{
#' \item{Disease: }{
#' \preformatted{
#' A disease defined in the database
#' OMIM, DECIPHER and/or Orphanet.}
#' }
#' \item{Phenotype: }{A clinical feature associated with one or more diseases.}
#' \item{Symptom: }{
#' \preformatted{
#' A phenotype within the context of a particular disease.
#' Within a given phenotype, there may be multiple symptoms with
#'  partially overlapping genetic mechanisms.}
#' }
#' \item{Assocation: }{
#' \preformatted{
#' A cell type-specific enrichment test result conducted
#' at the disease-level, phenotype-level, or symptom-level.}
#' }
#' }
#' @param keep_celltypes Cell type to keep.
#' @param sort_cols How to sort the rows using \link[data.table]{setorderv}.
#'  \code{names(sort_cols)} will be supplied to the \code{cols=} argument
#'  and values will be supplied to the \code{order=} argument.
#' @param top_n Top N genes to keep when grouping by \code{group_vars}.
#' @param group_vars Columns to group by when selecting \code{top_n} genes.
#' @param gene_size Min/max gene size (important for therapeutics design).
#' @param keep_biotypes Which gene biotypes to keep.
#'  (e.g. "protein_coding", "processed_transcript", "snRNA",
#'  "lincRNA", "snoRNA", "IG_C_gene")
#' @param keep_specificity_quantiles Which cell type
#' specificity quantiles to keep (max quantile is 40).
#' @param keep_mean_exp_quantiles Which cell type
#' mean expression quantiles to keep (max quantile is 40).
#' @param return_report If \code{TRUE}, will return a named list containing a
#' \code{report} that shows the number of
#'  phenotypes/celltypes/genes remaining after each filtering step.
#' @param symptom_intersection_threshold Minimum proportion of genes
#' overlapping between a symptom gene list
#' (phenotype-associated genes in the context of a particular disease)
#' and the phenotype-cell type association driver genes
#' @param severity_threshold_max The max severity score that a phenotype can
#' have across any disease.
#' @param severity_score_gpt_threshold The minimum GPT severity score that a
#'  phenotype can have across any disease.
#' @param run_prune_ancestors Prune redundant ancestral terms if any of their
#' descendants are present. Passes to \link[KGExplorer]{prune_ancestors}.
#' @param ctd_list A named list of CellTypeDataset objects each
#'  created with \link[EWCE]{generate_celltype_data}.
#' @param effect_var Name of the effect size column in the \code{results}.
#' @param info_content_threshold Minimum phenotype information content
#' threshold.
#' @inheritParams ewce_para
#' @inheritParams ggnetwork_plot_full
#' @inheritParams EWCE::bootstrap_enrichment_test
#' @inheritParams HPOExplorer::phenos_to_granges
#' @inheritParams HPOExplorer::add_
#' @returns A data.table of the prioritised phenotype- and
#'  cell type-specific gene targets.
#'
#' @export
#' @import HPOExplorer
#' @import data.table
#' @importFrom utils head
#' @examples
#' results = load_example_results()[seq(5000),]
#' out <- prioritise_targets(results=results)
prioritise_targets <- function(#### Input data ####
                               results = load_example_results(),
                               ctd_list = load_example_ctd(
                                 c("ctd_DescartesHuman.rds",
                                   "ctd_HumanCellLandscape.rds"),
                                 multi_dataset=TRUE),
                               phenotype_to_genes =
                                 HPOExplorer::load_phenotype_to_genes(),
                               hpo = HPOExplorer::get_hpo(),
                               #### Disease-level ####
                               keep_deaths =
                                 HPOExplorer::list_deaths(
                                   exclude=c("Miscarriage",
                                             "Stillbirth",
                                             "Prenatal death"),
                                   include_na = TRUE
                                 ),
                               #### Phenotype level ####
                               keep_descendants = c("Phenotypic abnormality"),
                               keep_ont_levels = NULL,
                               pheno_ndiseases_threshold = NULL,
                               gpt_filters = NULL,
                               severity_score_gpt_threshold=20,
                               keep_tiers = NULL,#c(1,2,NA),
                               severity_threshold_max = NULL,
                               info_content_threshold=8,
                               run_prune_ancestors=TRUE,
                               #### Symptom level ####
                               severity_threshold = NULL,#c(2,NA),
                               pheno_frequency_threshold = NULL,
                               keep_onsets =
                                 HPOExplorer::list_onsets(
                                   # exclude=c("Antenatal",
                                   #           "Fetal",
                                   #           "Congenital"),
                                   include_na = TRUE
                                 ),
                               #### Celltype level ####
                               effect_var="logFC",
                               q_threshold = 0.05,
                               effect_threshold = 1,
                               # symptom_p_threshold = NULL,
                               symptom_intersection_threshold = .25,
                               keep_celltypes = NULL,#terminal_celltypes()$CellType,
                               #### Gene level ####
                               evidence_score_threshold = 15,
                               keep_chr = c(seq(22),"X","Y"),
                               gene_size = list("min"=0,
                                                "max"=Inf),
                               gene_frequency_threshold = NULL,
                               keep_biotypes = NULL,
                               keep_specificity_quantiles = seq(30,40),
                               keep_mean_exp_quantiles = seq(30,40),
                               # symptom_gene_overlap = TRUE,
                               #### Sorting ####
                               sort_cols = c("severity_score_gpt"=-1,
                                             "q"=1,
                                             "logFC"=-1,
                                             "specificity"=-1,
                                             "mean_exp"=-1,
                                             "pheno_freq_mean"=-1,
                                             "gene_freq_mean"=-1,
                                             "width"=1),
                               top_n = NULL,
                               group_vars = c(
                                              # "disease_id",
                                              "hpo_id"
                                              # "CellType"
                                              ),
                               return_report = TRUE,
                               verbose = TRUE){
  q <- CellType <- width <- seqnames <- gene_biotype <-
    Severity_score <- cl_name <- cl_id <- Severity_score_max <-
    info_content <- severity_score_gpt <- NULL;

  force(results)
  force(ctd_list)
  force(hpo)
  force(phenotype_to_genes)
  t1 <- Sys.time()
  messager("Prioritising gene targets.",v=verbose)
  #### Add logFC ####
  add_logfc(results = results)
  #### add_hpo_id  #####
  results <- HPOExplorer::add_hpo_id(phenos = results,
                                     hpo = hpo)
  results <- HPOExplorer::add_hpo_name(phenos = results,
                                       hpo = hpo)
  #### add_hpo_definition  #####
  results <- HPOExplorer::add_hpo_definition(phenos = results,
                                             verbose = verbose)
  #### add_info_content #####
  results <- HPOExplorer::add_info_content(phenos = results,
                                           hpo = hpo)
  #### Add disease columns ####
  ## Add this early so I can get a count of the number of diseases
  ## in the final report plot.
  results <- HPOExplorer::add_disease(phenos = results,
                                      allow.cartesian = TRUE)
  #### start ####
  rep_dt <- report(dt = results,
                   step = "start",
                   verbose = verbose)

  #### Filter associations #####
  #### q_threshold ####
  if(!is.null(q_threshold)){
    messager("Filtering @ q-value <=",q_threshold,v=verbose)
    results <- results[q<=q_threshold,]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "q_threshold",
                   verbose = verbose)
  #### effect_threshold ####
  if(!is.null(effect_threshold)){
    messager("Filtering @",effect_var,">=",effect_threshold,v=verbose)
    results <- results[get(effect_var)>=effect_threshold,]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "effect_threshold",
                   verbose = verbose)
  # #### Filter symptoms ####
  # ## Do these steps early bc it will drastically reduce data size
  # ## and thus speed up all subsequent steps.
  # #### Filter smyptom overlap ####
  # results2 <- add_driver_genes(results = results,
  #                              ctd_list = ctd_list,
  #                              keep_quantiles = keep_specificity_quantiles)
  # #### symptom_p_threshold ####
  # if(!is.null(symptom_p_threshold)){
  #   results <- results[symptom.pval<symptom_p_threshold]
  # }
  # rep_dt <- report(dt = results,
  #                  rep_dt = rep_dt,
  #                  step = "symptom_p_threshold",
  #                  verbose = verbose)
  # #### symptom_intersection_threshold ####
  # if("intersection_size" %in% names(results)){
  #   if(!is.null(symptom_intersection_threshold)){
  #     results <- results[intersection_size>=symptom_intersection_threshold]
  #   }
  #   rep_dt <- report(dt = results,
  #                    rep_dt = rep_dt,
  #                    step = "symptom_intersection_threshold",
  #                    verbose = verbose)
  # }

  #### Filter diseases ####
  #### keep_deaths ####
  results <- HPOExplorer::add_death(phenos = results,
                                    keep_deaths = keep_deaths,
                                    agg_by = "disease_id",
                                    allow.cartesian = TRUE)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_deaths",
                   verbose = verbose)
  #### Filter phenotypes ####
  #### keep_descendants ####
  results <- HPOExplorer::add_ancestor(phenos = results,
                                       hpo = hpo,
                                       keep_descendants = keep_descendants)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_descendants",
                   verbose = verbose)
  #### keep_ont_levels ####
  results <- HPOExplorer::add_ont_lvl(phenos = results,
                                      absolute = TRUE,
                                      hpo = hpo,
                                      keep_ont_levels = keep_ont_levels)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_ont_levels",
                   verbose = verbose)
  #### gpt_filters ####
  ## Add GPT annotations
  results <- HPOExplorer::add_gpt_annotations(phenos = results,
                                              gpt_filters = gpt_filters)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "gpt_filters",
                   verbose = verbose)
  #### severity_score_gpt_threshold ####
  if(!is.null(severity_score_gpt_threshold)){
    results <- results[severity_score_gpt>=severity_score_gpt_threshold,]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "severity_score_gpt_threshold",
                   verbose = verbose)
  #### info_content_threshold ####
  if(!is.null(info_content_threshold)){
    results <- results[info_content>=info_content_threshold,]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "info_content_threshold",
                   verbose = verbose)

  #### keep_onsets ####
  results <- HPOExplorer::add_onset(phenos = results,
                                    keep_onsets = keep_onsets,
                                    agg_by=c("disease_id","hpo_id"))
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_onsets",
                   verbose = verbose)
  #### keep_tiers ####
  results <- HPOExplorer::add_tier(phenos = results,
                                   hpo = hpo,
                                   keep_tiers = keep_tiers,
                                   verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_tiers",
                   verbose = verbose)
  #### severity_threshold ####
  results <- HPOExplorer::add_severity(phenos = results,
                                       severity_threshold = severity_threshold)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "severity_threshold",
                   verbose = verbose)
  #### severity_threshold_max ####
  ## i.e. is a phenotype always severe, regardless of disease?
  if(!is.null(severity_threshold_max)){
    results <- results[,Severity_score_max:=gsub(
       -Inf,NA,max(Severity_score,na.rm = TRUE)),
       by="hpo_id"][Severity_score_max<=severity_threshold_max]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "severity_threshold_max",
                   verbose = verbose)
  #### pheno_ndiseases_threshold ####
  results <- HPOExplorer::add_ndisease(
    phenos = results,
    pheno_ndiseases_threshold = pheno_ndiseases_threshold,
    verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "pheno_ndiseases_threshold",
                   verbose = verbose)
  #### pheno_frequency_threshold ####
  results <- HPOExplorer::add_pheno_frequency(
    phenos = results,
    pheno_frequency_threshold = pheno_frequency_threshold)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "pheno_frequency_threshold",
                   verbose = verbose)

  #### Filter celltypes ####
  ### Fix celltypes
  results[,CellType:=EWCE::fix_celltype_names(CellType,
                                              make_unique = FALSE)]
  if(!is.null(keep_celltypes)){
    all_celltypes <- unique(results$CellType)
    results <- results[CellType %in% keep_celltypes|
                       cl_name %in% keep_celltypes|
                       cl_id %in% keep_celltypes,]
    valid_celltypes <- unique(results$CellType)
    messager(formatC(length(valid_celltypes),big.mark = ","),"/",
             formatC(length(all_celltypes)),
             "of cell types kept.",v=verbose)
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_celltypes",
                   verbose = verbose)
  #### Filter genes ####
  #### Add genes ####
  results <- HPOExplorer::add_genes(phenos=results,
                                    phenotype_to_genes=phenotype_to_genes,
                                    hpo = hpo)
  #### symptom_gene_overlap ####
  results <- HPOExplorer::phenos_to_granges(
    phenos = results,
    phenotype_to_genes = phenotype_to_genes,
    hpo = hpo,
    # gene_col = if(isTRUE(symptom_gene_overlap)) "intersection" else NULL,
    keep_chr = NULL,
    split.field = NULL,
    as_datatable = TRUE,
    verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "symptom_gene_overlap",
                   verbose = verbose)
  #### keep_chr ####
  if(!is.null(keep_chr)){
    messager("Filtering by keep_chr.",v=verbose)
    results <- results[seqnames %in% keep_chr,]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_chr",
                   verbose = verbose)
  #### evidence_score_threshold ####
  messager("Filtering by gene-disease association evidence.",
           v=verbose)
  results <- HPOExplorer::add_evidence(phenos = results,
                                       # allow.cartesian = TRUE,
                                       evidence_score_threshold = evidence_score_threshold)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "evidence_score_threshold",
                   verbose = verbose)
  #### gene_size ####
  if(!is.null(gene_size)){
    messager("Filtering by gene size.",v=verbose)
    ngenes <- length(unique(results$gene_symbol))
    results <- results[width>gene_size$min & width<gene_size$max,]
    messager(formatC(length(unique(results$gene_symbol)),big.mark = ","),"/",
             formatC(ngenes,big.mark = ","),"genes kept.",v=verbose)
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "gene_size",
                   verbose = verbose)
  #### keep_biotypes ####
  if(!is.null(keep_biotypes)){
    messager("Filtering by gene biotypes.",v=verbose)
    results <- results[gene_biotype %in% keep_biotypes,]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_biotypes",
                   verbose = verbose)
  ##### keep_specificity_quantiles ####
  ##### keep_mean_exp_quantiles ####
  results <- add_driver_genes(results = results,
                              ctd_list = ctd_list,
                              keep_quantiles = keep_specificity_quantiles,
                              metric = "specificity_quantiles",
                              group_var = c("hpo_id","disease_id")
                              )
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "add_driver_genes",
                   verbose = verbose)
  #### Symptom-level ####
  results <- add_symptom_results(results = results,
                                 q_threshold = q_threshold,
                                 celltype_col="CellType",
                                 ctd_list = ctd_list,
                                 phenotype_to_genes = phenotype_to_genes,
                                 keep_quantiles = keep_specificity_quantiles,
                                 top_n = NULL,
                                 drop_subthreshold=TRUE,
                                 proportion_driver_genes_symptom_threshold=symptom_intersection_threshold)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "symptom_intersection_threshold",
                   verbose = verbose)
  # ctd_out <- add_ctd(
  #   results = results,
  #   ctd = ctd,
  #   annotLevel = annotLevel,
  #   keep_specificity_quantiles = keep_specificity_quantiles,
  #   keep_mean_exp_quantiles = keep_mean_exp_quantiles,
  #   all.x = FALSE,
  #   rep_dt = rep_dt,
  #   verbose = verbose)
  # results <- ctd_out$results;
  # rep_dt <- ctd_out$rep_dt;
  #### gene_frequency_threshold ####
  results <- HPOExplorer::add_gene_frequency(
    phenotype_to_genes = results,
    gene_frequency_threshold = gene_frequency_threshold,
    allow.cartesian = TRUE,
    verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "gene_frequency_threshold",
                   verbose = verbose)

  #### prune_ancestors #####
  ## Even though this is a phenotype-level filter, do this at the end when we
  ## know which phenotype we have left.
  if(isTRUE(run_prune_ancestors)){
    results <- KGExplorer::prune_ancestors(dat = results,
                                           id_col = "hpo_id",
                                           ont = hpo)
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "prune_ancestors",
                   verbose = verbose)
  #### Sort genes ####
  # 1=ascending, -1=descending
  messager("Sorting rows.",v=verbose)
  sort_cols <- sort_cols[names(sort_cols) %in% names(results)]
  data.table::setorderv(results,
                        cols = names(sort_cols),
                        order = unname(sort_cols),
                        na.last = TRUE)
  #### top_n ####
  if(is.null(top_n)){
    top_targets <- results
  } else {
    messager("Finding top",top_n,"gene targets per:",
             paste(group_vars,collapse = ", "),v=verbose)
    top_targets <- results[,utils::head(.SD, top_n),
                             by = c(group_vars)]
    rep_dt <- report(dt = top_targets,
                     rep_dt = rep_dt,
                     step = "top_n",
                     verbose = verbose)
  }
  #### end ####
  rep_dt <- report(dt = top_targets,
                   rep_dt = rep_dt,
                   step = "end",
                   verbose = verbose)
  ## Add row diff
  rep_dt$Rows_diff <- c(0,
                        rep_dt$Rows[seq(nrow(rep_dt)-1)+1] -
                        rep_dt$Rows[seq(nrow(rep_dt)-1)])
  if(isTRUE(verbose)) round(difftime(Sys.time(),t1,units = "s"),0)
  #### Return ####
  if(isTRUE(return_report)){
    return(list(top_targets=top_targets,
                report=rep_dt))
  } else {
    return(top_targets)
  }
}

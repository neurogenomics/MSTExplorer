#' Prioritise target genes
#'
#' Prioritise target genes based on a procedure:\cr
#' \enumerate{
#' \item{Disease-level: \code{keep_deaths}: }{
#' Keep only diseases with a certain age of death.}
#' \item{Phenotype-level: \code{remove_descendants}: }{
#' Remove phenotypes belonging to a certain branch of the HPO,
#' as defined by an ancestor term.}
#' \item{Phenotype-level: \code{keep_ont_levels}: }{
#' Keep only phenotypes at certain absolute ontology levels within the HPO.}
#' \item{Phenotype-level: \code{pheno_ndiseases_threshold}: }{
#' The maximum number of diseases each phenotype can be associated with.}
#' \item{Phenotype-level: \code{keep_tiers}: }{
#' Keep only phenotypes with high severity Tiers.}
#' \item{Phenotype-level: \code{severity_threshold}: }{
#' Keep only phenotypes with mean Severity equal to or below the threshold.}
#' \item{Symptom-level: \code{pheno_frequency_threshold}: }{
#' Keep only phenotypes with mean frequency equal to or above the threshold
#'  (i.e. how frequently a phenotype is associated with any diseases in
#'  which it occurs).}
#' \item{Symptom-level: \code{keep_onsets}: }{
#' Keep only symptoms with a certain age of onset.}
#' \item{Symptom-level: \code{symptom_p_threshold}: }{
#' Uncorrected p-value threshold to filter cell type-symptom associations by.}
#' \item{Symptom-level: \code{symptom_intersection_size_threshold}: }{
#' Minimum number of genes overlapping between a symptom gene list
#' (phenotype-associated genes in the context of a particular disease)
#' and the celltype (genes in top 10/40 specificity quantiles).}
#' \item{Cell type-level: \code{q_threshold}: }{
#' Keep only cell type-phenotype association results at q<=0.05.}
#' \item{Cell type-level: \code{fold_threshold}: }{
#' Keep only cell type-phenotype association results at fold_change>=1.}
#' \item{Cell type-level: \code{keep_celltypes}: }{
#' Keep only terminally differentiated cell types.}
#' \item{Gene-level: \code{keep_seqnames}: }{
#' Remove genes on non-standard chromosomes.}
#' \item{Gene-level: \code{gene_size}: }{
#' Keep only genes <4.3kb in length.}
#' \item{Gene-level: \code{keep_biotypes}: }{
#' Keep only genes belonging to certain biotypes.}
#' \item{Gene-level: \code{gene_frequency_threshold}: }{
#' Keep only genes at or above a certain mean frequency threshold
#'  (i.e. how frequently a gene is associated with a given phenotype
#'  when observed within a disease).}
#' \item{Gene-level: \code{keep_specificity_quantiles}: }{
#'  Keep only genes in top specificity quantiles
#'  from the cell type dataset (\code{ctd}).}
#' \item{Gene-level: \code{keep_mean_exp_quantiles}: }{
#'  Keep only genes in top mean expression quantiles
#'  from the cell type dataset (\code{ctd}).}
#' \item{Gene-level: \code{symptom_gene_overlap}: }{
#'  Ensure that genes nominated at the phenotype-level also
#'  appear in the genes overlapping at the cell type-specific symptom-level.}
#' \item{All levels: \code{top_n}: }{
#' Sort candidate targets by a preferred order of metrics and
#'  only return the top N targets per cell type-phenotype combination.}
#' }
#'
#' Term key:\cr
#' \itemize{
#' \item{Disease: }{A disease defined in the database
#' OMIM, DECIPHER and/or Orphanet.}
#' \item{Phenotype: }{A clinical feature associated with one or more diseases.}
#' \item{Symptom: }{A phenotype within the context of a particular disease.
#' Within a given phenotype, there may be multiple symptoms with
#'  partially overlapping genetic mechanisms.}
#' \item{Assocation: }{A cell type-specific enrichment test result conducted
#' at the disease-level, phenotype-level, or symptom-level.}
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
#' @param symptom_p_threshold The p-value threshold of celltype-symptom
#' enrichment results (using  \link[MultiEWCE]{gen_overlap}).
#' Here, "symptoms" are defined as the presentation of a phenotype in the
#' context of a particular disease.
#'  In other words:
#'  phenotype (HPO_ID) + disease (DatabaseID) = symptom (HPO_ID.DatabaseID)
#' @param symptom_intersection_size_threshold
#' The minimum number of intersecting genes between a symptom and a celltype to
#' consider it a significant enrichment. Refers to the result from
#' \link[MultiEWCE]{gen_overlap}.
#' @param symptom_gene_overlap The gene for a particular symptom
#' (phenotype + disease) must appear in the celltype-symptom
#'  enrichment results.
#' @inheritParams ewce_para
#' @inheritParams ggnetwork_plot_full
#' @inheritParams EWCE::bootstrap_enrichment_test
#' @inheritParams HPOExplorer::phenos_to_granges
#' @inheritParams HPOExplorer::add_ancestor
#' @inheritParams HPOExplorer::add_severity
#' @inheritParams HPOExplorer::add_tier
#' @inheritParams HPOExplorer::add_death
#' @inheritParams HPOExplorer::add_onset
#' @inheritParams HPOExplorer::add_gene_frequency
#' @inheritParams HPOExplorer::add_pheno_frequency
#' @inheritParams HPOExplorer::add_ndisease
#' @inheritParams HPOExplorer::add_ont_lvl
#' @returns A data.table of the prioritised phenotype- and
#'  celltype-specific gene targets.
#'
#' @export
#' @import HPOExplorer
#' @import data.table
#' @importFrom utils head
#' @examples
#' res <- prioritise_targets(symptom_p_threshold = NULL)
prioritise_targets <- function(#### Input data ####
                               results = load_example_results(),
                               ctd = load_example_ctd(),
                               annotLevel = 1,
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
                               remove_descendants = c("Clinical course"),
                               keep_ont_levels = NULL,
                               pheno_ndiseases_threshold = NULL,
                               keep_tiers = c(1,2,NA),
                               severity_threshold = c(2,NA),
                               #### Symptom level ####
                               pheno_frequency_threshold = NULL,
                               keep_onsets =
                                 HPOExplorer::list_onsets(
                                   exclude=c("Antenatal",
                                             "Fetal",
                                             "Congenital"),
                                   include_na = TRUE
                                 ),
                               symptom_p_threshold = NULL,
                               symptom_intersection_size_threshold = 1,
                               #### Celltype level ####
                               q_threshold = 0.05,
                               fold_threshold = 1,
                               keep_celltypes = terminal_celltypes()$CellType,
                               #### Gene level ####
                               keep_seqnames = c(seq_len(22),"X","Y"),
                               gene_size = list("min"=0,
                                                "max"=4300),
                               gene_frequency_threshold = NULL,
                               keep_biotypes = NULL,
                               keep_specificity_quantiles = NULL,
                               keep_mean_exp_quantiles = seq(1,40),
                               symptom_gene_overlap = TRUE,
                               #### Sorting ####
                               sort_cols = c("tier_merge"=1,
                                             "Severity_score_mean"=1,
                                             "q"=1,
                                             "fold_change"=-1,
                                             "specificity"=-1,
                                             "mean_exp"=-1,
                                             "pheno_freq_mean"=-1,
                                             "gene_freq_mean"=-1,
                                             "width"=1),
                               top_n = NULL,
                               group_vars = c("HPO_ID","CellType"),
                               return_report = TRUE,
                               verbose = TRUE){

  # o <- devoptera::args2vars(prioritise_targets)

  q <- fold_change <- CellType <- Gene <- HPO_ID <-
    HPO_term_valid <- specificity_quantile <- intersection <- symptom.pval <-
    mean_exp_quantile <- celltype_fixed <- intersection_size <- NULL;

  t1 <- Sys.time()
  messager("Prioritising gene targets.",v=verbose)
  #### add_hpo_id  #####
  results <- HPOExplorer::add_hpo_id(phenos = results,
                                     phenotype_to_genes = phenotype_to_genes,
                                     hpo = hpo,
                                     verbose = verbose)
  if("HPO_term_valid" %in% names(results)){
    results <- results[HPO_term_valid==TRUE,]
  }
  #### add_hpo_definition  #####
  results <- HPOExplorer::add_hpo_definition(phenos = results,
                                             verbose = verbose)
  #### add_info_content #####
  if("info_content" %in% names(sort_cols)){
    results <- HPOExplorer::add_info_content(phenos = results,
                                             hpo = hpo,
                                             verbose = verbose)
  }
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
  #### fold_threshold ####
  if(!is.null(fold_threshold)){
    messager("Filtering @ fold-change >=",fold_threshold,v=verbose)
    results <- results[fold_change>=fold_threshold,]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "fold_threshold",
                   verbose = verbose)
  #### Filter phenotypes ####
  #### remove_descendants ####
  results <- HPOExplorer::add_ancestor(phenos = results,
                                       hpo = hpo,
                                       remove_descendants = remove_descendants,
                                       verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "remove_descendants",
                   verbose = verbose)
  #### keep_ont_levels ####
  results <- HPOExplorer::add_ont_lvl(phenos = results,
                                      absolute = TRUE,
                                      keep_ont_levels = keep_ont_levels,
                                      verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_ont_levels",
                   verbose = verbose)
  #### keep_onsets ####
  results <- HPOExplorer::add_onset(phenos = results,
                                    keep_onsets = keep_onsets,
                                    agg_by=c("DatabaseID","HPO_ID"),
                                    verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_onsets",
                   verbose = verbose)
  #### keep_deaths ####
  results <- HPOExplorer::add_death(phenos = results,
                                    keep_deaths = keep_deaths,
                                    agg_by = "DatabaseID",
                                    verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_deaths",
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
                                       severity_threshold = severity_threshold,
                                       verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "severity_threshold",
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
    pheno_frequency_threshold = pheno_frequency_threshold,
    verbose = verbose)
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "pheno_frequency_threshold",
                   verbose = verbose)
  #### symptom_p_threshold ####
  if(!is.null(symptom_p_threshold)){
    results <- results[symptom.pval<symptom_p_threshold]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "symptom_p_threshold",
                   verbose = verbose)
  #### symptom_intersection_size_threshold ####
  if(!is.null(symptom_intersection_size_threshold)){
    results <- results[intersection_size>=symptom_intersection_size_threshold]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "symptom_intersection_size_threshold",
                   verbose = verbose)


  #### Filter celltypes ####
  if(!is.null(keep_celltypes)){
    all_celltypes <- unique(results$CellType)
    results <- results[CellType %in% keep_celltypes,]
    valid_celltypes <- unique(results$CellType)
    messager(formatC(length(valid_celltypes),big.mark = ","),"/",
             formatC(length(all_celltypes)),
             "of cell types kept.",v=verbose)
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_celltypes",
                   verbose = verbose)
  #### Filter genes by size ####
  messager("Filtering by gene size.",v=verbose)
  gr <- HPOExplorer::phenos_to_granges(phenos = results,
                                       phenotype_to_genes = phenotype_to_genes,
                                       hpo = hpo,
                                       keep_seqnames = keep_seqnames,
                                       split.field = NULL,
                                       verbose = verbose)
  rep_dt <- report(dt = gr,
                   rep_dt = rep_dt,
                   step = "keep_seqnames",
                   verbose = verbose)
  #### gene_size ####
  if(!is.null(gene_size)){
    ngenes <- length(unique(gr$Gene))
    gr <- gr[gr@ranges@width>gene_size$min & gr@ranges@width<gene_size$max,]
    messager(formatC(length(unique(gr$Gene)),big.mark = ","),"/",
             formatC(ngenes,big.mark = ","),"genes kept.",v=verbose)
  }
  rep_dt <- report(dt = gr,
                   rep_dt = rep_dt,
                   step = "gene_size",
                   verbose = verbose)
  #### keep_biotypese ####
  if(!is.null(keep_biotypes)){
    messager("Filtering by gene biotypes.",v=verbose)
    gr <- gr[gr$gene_biotype %in% keep_biotypes,]
  }
  rep_dt <- report(dt = gr,
                   rep_dt = rep_dt,
                   step = "keep_biotypes",
                   verbose = verbose)
  #### Identify high-specificity genes ####
  shared_genes <- intersect(gr$Gene,rownames(ctd[[annotLevel]]$specificity))
  specq_df <- make_specificity_dt(ctd = ctd,
                                 annotLevel = annotLevel,
                                 shared_genes = shared_genes,
                                 metric = "specificity_quantiles")
  spec_df <- make_specificity_dt(ctd = ctd,
                                 annotLevel = annotLevel,
                                 shared_genes = shared_genes,
                                 metric = "specificity")
  exp_df <- make_specificity_dt(ctd = ctd,
                                 annotLevel = annotLevel,
                                 shared_genes = shared_genes,
                                 metric = "mean_exp")
  expq_df <- make_specificity_dt(ctd = ctd,
                                annotLevel = annotLevel,
                                shared_genes = shared_genes,
                                metric = "mean_exp_quantiles")
  ##### Filter by genes  previously identified ####
  specq_df <- specq_df[Gene %in% gr$Gene,]
  expq_df <- expq_df[Gene %in% gr$Gene,]
  shared_genes <- intersect(intersect(specq_df$Gene,
                                      expq_df$Gene),
                            shared_genes)
  #### Filter by specificity #####
  if(!is.null(keep_specificity_quantiles)){
    messager("Filtering by specificity_quantile.",v=verbose)
    specq_df <- specq_df[specificity_quantile %in% keep_specificity_quantiles,]
  }
  #### Filter by specificity #####
  if(!is.null(keep_mean_exp_quantiles)){
    messager("Filtering by mean_exp_quantile.",v=verbose)
    expq_df <- expq_df[mean_exp_quantile %in% keep_mean_exp_quantiles,]
  }
  ##### Filter by cell types  previously identified ####
  keep_celltypes2 <- unique(results$CellType)[
    EWCE::fix_celltype_names(unique(results$CellType)) %in%
      specq_df$celltype_fixed
  ]
  results <- results[CellType %in% keep_celltypes2,]
  #### Merge genes with phenotype/celltype results ####
  gr_df <- data.table::as.data.table(gr)[Gene %in% shared_genes,]
  data.table::setkeyv(gr_df,"HPO_ID")
  cols <- unique(
    c("Phenotype","HPO_ID",
      names(results)[!names(results) %in% names(gr_df)])
  )
  df_merged <- data.table::merge.data.table(
    x = results[HPO_ID %in% unique(gr_df$HPO_ID),][,cols,with=FALSE],
    y = unique(
      gr_df[,c("HPO_ID","Gene","gene_biotype",
               "seqnames","start","end","width")]
    ),
    all = TRUE,
    allow.cartesian = TRUE,
    by = "HPO_ID"
  )
  ct_dict <- stats::setNames(
    EWCE::fix_celltype_names(celltypes = unique(df_merged$CellType)),
    unique(df_merged$CellType))
  df_merged[,celltype_fixed:=ct_dict[CellType]]
  by_cols <- c("Gene","celltype_fixed")
  #### gene_frequency_threshold ####
  df_merged <- HPOExplorer::add_gene_frequency(
    phenotype_to_genes = df_merged,
    gene_frequency_threshold = gene_frequency_threshold,
    verbose = verbose)
  rep_dt <- report(dt = df_merged,
                   rep_dt = rep_dt,
                   step = "gene_frequency_threshold",
                   verbose = verbose)
  #### Merge: specificity ####
  df_merged <- df_merged |>
    data.table::merge.data.table(y = spec_df,
                                 by = by_cols)
  #### Merge: specificity_quantiles ####
  df_merged <- df_merged |>
    data.table::merge.data.table(y = specq_df,
                                 by = by_cols)
  rep_dt <- report(dt = df_merged,
                   rep_dt = rep_dt,
                   step = "keep_specificity_quantiles",
                   verbose = verbose)
  #### Merge: mean_exp ####
  df_merged <- df_merged |>
    data.table::merge.data.table(y = exp_df,
                                 by = by_cols)
  #### Merge: mean_exp_quantiles ####
  df_merged <- df_merged |>
    data.table::merge.data.table(y = expq_df,
                                 by = by_cols)
  rep_dt <- report(dt = df_merged,
                   rep_dt = rep_dt,
                   step = "keep_mean_exp_quantiles",
                   verbose = verbose)
  #### Filter by symptom intersection genes ####
  if(isTRUE(symptom_gene_overlap) &&
     "intersection" %in% names(df_merged)){
    df_merged <-
      df_merged[,symptom.in_intersection:=Gene %in% unlist(intersection),
                by="HPO_ID.LinkID"][symptom.in_intersection==TRUE,]
  }
  rep_dt <- report(dt = df_merged,
                   rep_dt = rep_dt,
                   step = "symptom_gene_overlap",
                   verbose = verbose)
  #### Sort genes ####
  # 1=ascending, -1=descending
  messager("Sorting rows.",v=verbose)
  sort_cols <- sort_cols[names(sort_cols) %in% names(df_merged)]
  data.table::setorderv(df_merged,
                        cols = names(sort_cols),
                        order = unname(sort_cols),
                        na.last = TRUE)

  if(is.null(top_n)){
    top_targets <- df_merged
  } else {
    #### top_n ####
    messager("Finding top",top_n,"gene targets per:",
             paste(group_vars,collapse = ", "),v=verbose)
    top_targets <- df_merged[,utils::head(.SD, top_n),
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
  if(isTRUE(verbose)) round(difftime(Sys.time(),t1,units = "s"),0)
  #### Return ####
  if(isTRUE(return_report)){
    return(list(top_targets=top_targets,
                report=rep_dt))
  } else {
    return(top_targets)
  }
}

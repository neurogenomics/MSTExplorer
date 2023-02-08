#' Prioritise target genes
#'
#' Prioritise target genes based on a procedure:
#' \enumerate{
#' \item{Association-level: \code{q_threshold}: }{
#' Keep only results at q<=0.05.}
#' \item{Association-level: \code{fold_threshold}: }{
#' Keep only results at fold_change>=1.}
#' \item{Phenotype-level: \code{keep_ont_levels}: }{
#' Keep only phenotypes at certain absolute ontology levels within the HPO.}
#' \item{Phenotype-level: \code{keep_onsets}: }{
#' Keep only phenotypes with certain age of onsets.}
#' \item{Phenotype-level: \code{keep_tiers}: }{
#' Keep only phenotypes with high severity Tiers.}
#' \item{Phenotype-level: \code{severity_threshold}: }{
#' Keep only phenotypes with mean Severity equal to or below the threshold.}
#' \item{Phenotype-level: \code{pheno_frequency_threshold}: }{
#' Keep only phenotypes with mean frequency equal to or above the threshold
#'  (i.e. how frequently a phenotype is associated with any diseases in
#'  which it occurs).}
#' \item{Cell type-level:  \code{keep_celltypes}: }{
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
#' \item{All levels: \code{top_n}: }{
#' Sort candidate targets by a preferred order of metrics and
#'  only return the top N targets per cell type-phenotype combination.}
#' }
#'
#' @param keep_celltypes Cell type to keep.
#' @param keep_tiers Tiers from \link[HPOExplorer]{hpo_tiers} to keep.
#'  Include \code{NA} if you wish to retain phenotypes that
#'  do not have any Tier assignment.
#' @param severity_threshold Only keep phenotypes with severity scores below the
#'  set threshold. The severity score ranges from 1-4 where 1 is the MOST severe.
#'  Include \code{NA} if you wish to retain phenotypes that
#'  do not have any severity score.
#' @param keep_ont_levels Only keep phenotypes at certain \emph{absolute}
#'  ontology levels to keep.
#' See \link[HPOExplorer]{add_ont_lvl} for details.
#' @param keep_onsets The age of onset associated with each HPO ID to keep.
#'  If >1 age of onset is associated with the term,
#'  only the earliest onset is considered.
#'  See \link[HPOExplorer]{add_onset} for details.
#' @param pheno_frequency_threshold Only keep phenotypes with frequency
#'  above the set threshold. Frequency ranges from 0-100 where 100 is
#'  a phenotype that occurs 100% of the time in all associated diseases.
#'  Include \code{NA} if you wish to retain phenotypes that
#'  do not have any frequency data.
#' @param gene_frequency_threshold Only keep genes with frequency
#'  above the set threshold. Frequency ranges from 0-100 where 100 is
#'  a gene that occurs 100% of the time in a given phenotype.
#'  Include \code{NA} if you wish to retain genes that
#'  do not have any frequency data.
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
#'
#' @inheritParams ewce_para
#' @inheritParams ggnetwork_plot_full
#' @inheritParams EWCE::bootstrap_enrichment_test
#' @inheritParams HPOExplorer::phenos_to_granges
#' @returns A data.table of the prioritised phenotype- and
#'  celltype-specific gene targets.
#'
#' @export
#' @importFrom HPOExplorer load_phenotype_to_genes get_hpo add_hpo_id
#' @importFrom HPOExplorer list_onsets phenos_to_granges add_tier add_onset
#' @importFrom HPOExplorer add_info_content add_ancestor
#' @importFrom data.table .SD := merge.data.table setorderv as.data.table
#' @importFrom data.table data.table
#' @importFrom utils head
#' @examples
#' res <- prioritise_targets()
prioritise_targets <- function(results = load_example_results(),
                               ctd = load_example_ctd(),
                               annotLevel = 1,
                               q_threshold = 0.05,
                               fold_threshold = 1,
                               keep_ont_levels = NULL,
                               keep_onsets =
                                 HPOExplorer::list_onsets(
                                   exclude_onsets=c("Antenatal","Fetal")
                                 ),
                               keep_tiers = c(1,2),
                               severity_threshold = c(2,NA),
                               pheno_frequency_threshold = 25,
                               keep_celltypes = terminal_celltypes()$CellType,
                               keep_seqnames = c(seq_len(22),"X","Y"),
                               gene_size = list("min"=0,
                                                "max"=4300),
                               gene_frequency_threshold = c(10,NA),
                               keep_biotypes = NULL,
                               keep_specificity_quantiles =
                                 seq(39,40),
                               keep_mean_exp_quantiles =
                                 keep_specificity_quantiles,
                               sort_cols = c("tier"=1,
                                             "tier_auto"=1,
                                             "Severity_score_mean"=1,
                                             "q"=1,
                                             "fold_change"=-1,
                                             "specificity_quantile"=-1,
                                             "mean_exp_quantile"=-1,
                                             "specificity"=-1,
                                             "mean_exp"=-1,
                                             "pheno_freq_mean"=-1,
                                             "gene_freq_mean"=-1,
                                             "width"=1),
                               top_n = 20,
                               group_vars = c("HPO_ID","CellType"),
                               phenotype_to_genes =
                                         HPOExplorer::load_phenotype_to_genes(),
                               hpo = HPOExplorer::get_hpo(),
                               return_report = TRUE,
                               verbose = TRUE){

  # templateR:::args2vars(prioritise_targets)

  q <- fold_change <- CellType <- Gene <- tier <- HPO_ID <-
    HPO_term_valid <- Onset_earliest <- specificity_quantile <-
    mean_exp_quantile <- celltype_fixed <- Severity_score_min <-
    tier_auto <- pheno_freq_mean <- gene_freq_mean <- ontLvl <- NULL;

  t1 <- Sys.time()
  messager("Prioritising gene targets.",v=verbose)
  #### add_hpo_id  #####
  results <- HPOExplorer::add_hpo_id(phenos = results,
                                     phenotype_to_genes = phenotype_to_genes,
                                     hpo = hpo)
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
  #### add_ancestor ####
  results <- HPOExplorer::add_ancestor(phenos = results,
                                       hpo = hpo,
                                       verbose = verbose)
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
  #### keep_ont_levels ####
  results <- HPOExplorer::add_ont_lvl(phenos = results,
                                      absolute = TRUE,
                                      verbose = verbose)
  if(!is.null(keep_ont_levels)){
    results <- results[ontLvl %in% keep_ont_levels,]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_ont_levels",
                   verbose = verbose)
  #### keep_onsets ####
  results <- HPOExplorer::add_onset(phenos = results,
                                    verbose = verbose)
  if(!is.null(keep_onsets)){
    results <- results[Onset_earliest %in% keep_onsets,]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_onsets",
                   verbose = verbose)
  #### keep_tiers ####
  results <- HPOExplorer::add_tier(phenos = results,
                                   hpo = hpo,
                                   verbose = verbose)
  if(!is.null(keep_tiers)){
    results <- results[(tier %in% keep_tiers) |
                       (tier_auto %in% keep_tiers),]
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "keep_tiers",
                   verbose = verbose)
  #### severity_threshold ####
  results <- HPOExplorer::add_modifier(phenos = results,
                                       verbose = verbose)
  if(!is.null(severity_threshold)){
    if(any(is.na(severity_threshold))){
      results <- results[
        Severity_score_min<=min(severity_threshold,na.rm = TRUE) |
        is.na(Severity_score_min),]
    } else{
      results <- results[Severity_score_min<=
                           min(severity_threshold,na.rm = TRUE),]
    }
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "severity_threshold",
                   verbose = verbose)
  #### pheno_frequency_threshold ####
  results <- HPOExplorer::add_pheno_frequency(phenos = results,
                                              verbose = verbose)
  if(!is.null(pheno_frequency_threshold)){
    if(any(is.na(pheno_frequency_threshold))){
      results <- results[
        pheno_freq_mean>=min(pheno_frequency_threshold,na.rm = TRUE) |
          is.na(pheno_freq_mean),]
    } else{
      results <- results[pheno_freq_mean>=
                           min(pheno_frequency_threshold,na.rm = TRUE),]
    }
  }
  rep_dt <- report(dt = results,
                   rep_dt = rep_dt,
                   step = "pheno_frequency_threshold",
                   verbose = verbose)
  #### Filter celltypes ####
  if(!is.null(keep_celltypes)){
    all_celltypes <- unique(results$CellType)
    results <- results[tolower(CellType) %in% tolower(keep_celltypes),]
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
  data.table::setnames(gr_df,"ID","HPO_ID")
  data.table::setkeyv(gr_df,"HPO_ID")
  cols <- unique(
    c("Phenotype","HPO_ID",
      names(results)[!names(results) %in% names(gr_df)])
  )
  df_merged <- data.table::merge.data.table(
    x = unique(
      results[HPO_ID %in% unique(gr_df$HPO_ID),][,cols,with=FALSE]
    ),
    y = unique(
      gr_df[,c("HPO_ID","Gene","gene_biotype",
               "seqnames","start","end","width")]
    ),
    all = TRUE,
    allow.cartesian = TRUE,
    by="HPO_ID"
  )
  ct_dict <- stats::setNames(
    EWCE::fix_celltype_names(celltypes = unique(df_merged$CellType)),
    unique(df_merged$CellType))
  df_merged[,celltype_fixed:=ct_dict[CellType]]
  by_cols <- c("Gene","celltype_fixed")
  #### gene_frequency_threshold ####
  df_merged <- HPOExplorer::add_gene_frequency(phenotype_to_genes = df_merged,
                                               verbose = verbose)
  if(!is.null(gene_frequency_threshold)){
    if(any(is.na(gene_frequency_threshold))){
      df_merged <- df_merged[
        gene_freq_mean>=min(gene_frequency_threshold,na.rm = TRUE) |
          is.na(gene_freq_mean),]
    } else{
      df_merged <- df_merged[gene_freq_mean>=
                               min(gene_frequency_threshold,na.rm = TRUE),]
    }
  }
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

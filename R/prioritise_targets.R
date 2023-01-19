#' Prioritise target genes
#'
#' Prioritise target genes based on a procedure:
#' \enumerate{
#' \item{Association-level: }{Keep only results at q<=0.05.}
#' \item{Association-level: }{Keep only results at fold_change>=1.}
#' \item{Phenotype-level: }{Keep only phenotypes with high severity Tiers}
#' \item{Phenotype-level: }{Keep only phenotypes with postnatal age of onsets.}
#' \item{Cell type-level: }{Keep only terminally differentiated cell types.}
#' \item{Gene-level: }{Remove genes on non-standard chromosomes.}
#' \item{Gene-level: }{Keep only genes <4.3kb in length.}
#' \item{Gene-level: }{Keep only genes in top specificity quantiles
#' from the cell type dataset (\code{ctd}).}
#' \item{Gene-level: }{Sort genes by cell type specificity and mean expression
#' from the cell type dataset (\code{ctd}).}
#' }
#' @param keep_celltypes Cell type to keep.
#' @param keep_tiers Tiers from \link[HPOExplorer]{hpo_tiers} to keep.
#' @param keep_onsets The age of onset associated with each HPO ID to keep.
#' @param sort_cols How to sort the rows using \link[data.table]{setorderv}.
#' \code{names(sort_cols)} will be supplied to the \code{cols=} argument
#' and values will be supplied to the \code{order=} argument.
#' @param top_n Top N genes to keep when grouping by \code{group_vars}.
#' @param group_vars Columns to group by when selecting \code{top_n} genes.
#' @param gene_size Min/max gene size (important for therapeutics design).
#' @param keep_biotypes Which gene biotypes to keep.
#' @param keep_specificity_quantiles Which cell type
#' specificity quantiles to keep (max specificity quantile is 40).
#' @inheritParams ewce_para
#' @inheritParams ggnetwork_plot_full
#' @inheritParams EWCE::bootstrap_enrichment_test
#' @inheritParams HPOExplorer::phenos_to_granges
#'
#' @export
#' @importFrom HPOExplorer load_phenotype_to_genes get_hpo add_hpo_id
#' @importFrom HPOExplorer list_onsets phenos_to_granges add_tier add_onset
#' @importFrom HPOExplorer add_info_content
#' @importFrom data.table .SD := merge.data.table setorderv as.data.table
#' @importFrom data.table data.table
#' @importFrom utils head
#' @examples
#' results <- load_example_results()
#' ctd <- load_example_ctd()
#' top_targets <- prioritise_targets(results = results,
#'                                   ctd = ctd)
prioritise_targets <- function(results,
                               ctd,
                               annotLevel = 1,
                               q_threshold = 0.05,
                               fold_threshold = 1,
                               keep_celltypes = terminal_celltypes()$CellType,
                               keep_onsets =
                                 HPOExplorer::list_onsets(postnatal_only=TRUE),
                               keep_tiers = 1,
                               gene_size = list("min"=0,
                                                "max"=4300),
                               keep_seqnames = c(seq_len(22),"X","Y"),
                               keep_biotypes = NULL,
                               keep_specificity_quantiles = seq(38,40),
                               sort_cols = c("tier"=1,
                                             "q"=1,
                                             "fold_change"=-1,
                                             "specificity_quantile"=-1,
                                             "specificity"=-1,
                                             "mean_exp"=-1,
                                             "width"=1),
                               top_n = 3,
                               group_vars = c("HPO_ID","CellType"),
                               phenotype_to_genes =
                                         HPOExplorer::load_phenotype_to_genes(),
                               hpo = HPOExplorer::get_hpo(),
                               verbose = TRUE){
  # templateR:::source_all()
  # templateR:::args2vars(prioritise_targets)

  q <- fold_change <- CellType <- Gene <- tier <- HPO_ID <-
    HPO_term_valid <- Onset <- specificity_quantile <- celltype_fixed <- NULL;

  t1 <- Sys.time()
  messager("Prioritising gene targets.",v=verbose)
  results <- HPOExplorer::add_hpo_id(phenos = results,
                                     phenotype_to_genes = phenotype_to_genes,
                                     hpo = hpo)
  results <- results[HPO_term_valid==TRUE,]
  #### Add info content #####
  if("info_content" %in% names(sort_cols)){
    results <- HPOExplorer::add_info_content(phenos = results,
                                             hpo = hpo,
                                             verbose = verbose)
  }
  report(dt = results,
         verbose = verbose)
  #### Filter associations #####
  if(!is.null(q_threshold)){
    messager("Filtering @ q-value <=",q_threshold,v=verbose)
    results <- results[q<=q_threshold,]
    report(dt = results,
           verbose = verbose)
  }
  if(!is.null(fold_threshold)){
    messager("Filtering @ fold-change >=",fold_threshold,v=verbose)
    results <- results[fold_change>=fold_threshold,]
    report(dt = results,
           verbose = verbose)
  }
  #### Filter phenotypes ####
  if(!is.null(keep_tiers)){
    results <- HPOExplorer::add_tier(phenos = results,
                                     all.x = FALSE,
                                     verbose = verbose)
    results <- results[tier %in% keep_tiers,]
    report(dt = results,
           verbose = verbose)
  }
  if(!is.null(keep_onsets)){
    results <- HPOExplorer::add_onset(phenos = results,
                                      all.x = TRUE,
                                      verbose = verbose)
    results <- results[Onset %in% keep_onsets,]
    report(dt = results,
           verbose = verbose)
  }
  #### Filter celltypes ####
  if(!is.null(keep_celltypes)){
    all_celltypes <- unique(results$CellType)
    results <- results[tolower(CellType) %in% tolower(keep_celltypes),]
    valid_celltypes <- unique(results$CellType)
    messager(formatC(length(valid_celltypes),big.mark = ","),"/",
             formatC(length(all_celltypes)),
             "of cell types kept.",v=verbose)
    report(dt = results,
           verbose = verbose)
  }
  #### Filter genes by size ####
  messager("Filtering by gene size.",v=verbose)
  gr <- HPOExplorer::phenos_to_granges(phenos = results,
                                       phenotype_to_genes = phenotype_to_genes,
                                       hpo = hpo,
                                       keep_seqnames = keep_seqnames,
                                       split.field = NULL,
                                       verbose = verbose)
  ngenes <- length(unique(gr$Gene))
  gr <- gr[gr@ranges@width>gene_size$min & gr@ranges@width<gene_size$max,]
  messager(formatC(length(unique(gr$Gene)),big.mark = ","),"/",
           formatC(ngenes,big.mark = ","),"genes kept.",v=verbose)
  #### Filter genes by biotype ####
  if(!is.null(keep_biotypes)){
    messager("Filtering by gene biotypes.",v=verbose)
    gr <- gr[gr$gene_biotype %in% keep_biotypes,]
  }
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

  ##### Filter by genes  previously identified ####
  specq_df <- specq_df[Gene %in% gr$Gene,]
  ##### Filter by cell types  previously identified ####
  keep_celltypes2 <- unique(results$CellType)[
    EWCE::fix_celltype_names(unique(results$CellType)) %in%
      specq_df$celltype_fixed
  ]
  results <- results[CellType %in% keep_celltypes2,]
  #### Filter by specificity #####
  messager("Filtering by specificity_quantile.",v=verbose)
  if(!is.null(keep_specificity_quantiles)){
    specq_df <- specq_df[specificity_quantile %in% keep_specificity_quantiles,]
  }
  #### Merge with main results ####
  gr_df <- data.table::as.data.table(gr)[Gene %in% shared_genes,]
  data.table::setnames(gr_df,"ID","HPO_ID")
  data.table::setkeyv(gr_df,"HPO_ID")
  report(dt = gr_df,
         verbose = verbose)
  cols <- unique(
    c("Phenotype","HPO_ID",
      names(results)[!names(results) %in% names(gr_df)])
  )
  df_merged <- data.table::merge.data.table(
    x = unique(
      results[HPO_ID %in% unique(gr_df$HPO_ID),][,cols,with=FALSE]
    ),
    y = unique(
      gr_df[,c("HPO_ID","Gene","gene_biotype","seqnames","start","end","width")]
    ),
    all = TRUE,
    allow.cartesian = TRUE,
    by="HPO_ID"
  )
  ct_dict <- stats::setNames(
    EWCE::fix_celltype_names(celltypes = unique(df_merged$CellType)),
    unique(df_merged$CellType))
  df_merged[,celltype_fixed:=ct_dict[CellType]]
  #### Add specificity / mean expression metrics ####
  df_merged <- df_merged |>
    data.table::merge.data.table(y = specq_df,
                                 by = c("Gene","celltype_fixed")) |>
    data.table::merge.data.table(y = spec_df,
                                 by = c("Gene","celltype_fixed")) |>
    data.table::merge.data.table(y = exp_df,
                                 by = c("Gene","celltype_fixed"))
  report(dt = df_merged,
         verbose = verbose)
  #### Sort genes ####
  # 1=ascending, -1=descending
  messager("Sorting rows.",v=verbose)
  data.table::setorderv(df_merged,
                        cols = names(sort_cols),
                        order = unname(sort_cols),
                        na.last = TRUE)

  if(is.null(top_n)){
    top_targets <- df_merged
  } else {
    messager("Finding top",top_n,"gene targets per:",
             paste(group_vars,collapse = ", "),v=verbose)
    top_targets <- df_merged[,utils::head(.SD, top_n),
                             by = c(group_vars)]
  }
  #### Report ####
  report(dt = top_targets,
         verbose = verbose)
  if(isTRUE(verbose)) round(difftime(Sys.time(),t1,units = "s"),0)
  #### Return ####
  return(top_targets)
}

#' Prioritise target genes
#'
#' Prioritise target genes based on a procedure:
#' \enumerate{
#' \item{Include only results at q<=0.05.}
#' \item{Include only results at fold_change>=1.}
#' \item{Include only terminally differentiated cells.}
#' }
#' @param keep_celltypes Cell type to keep.
#' @param keep_tiers Tiers from \link[HPOExplorer]{hpo_tiers} to keep.
#' @inheritParams ewce_para
#' @inheritParams ggnetwork_plot_full
#' @inheritParams EWCE::bootstrap_enrichment_test
#'
#' @export
#' @importFrom HPOExplorer load_phenotype_to_genes get_hpo add_hpo_id
#' @examples
#' results <- load_example_results()
#' ctd <- load_example_ctd()
priortise_targets <- function(results,
                              ctd,
                              annotLevel = 1,
                              q_threshold = .05,
                              fold_threshold = 1,
                              keep_celltypes = NULL,
                              keep_tiers = c(1,2),
                              phenotype_to_genes =
                                        HPOExplorer::load_phenotype_to_genes(),
                              hpo = HPOExplorer::get_hpo()){
  # templateR:::source_all()
  # templateR:::args2vars(priortise_target)

  q <- fold_change <- CellType <- Phenotype <- Gene <- tier <- HPO_ID <- NULL;

  results <- HPOExplorer::add_hpo_id(phenos = results,
                                     phenotype_to_genes = phenotype_to_genes,
                                     hpo = hpo)
  #### Filter associations #####
  if(!is.null(q_threshold)){
    results <- results[q<=q_threshold,]
  }
  if(!is.null(fold_threshold)){
    results <- results[fold_change>=fold_threshold,]
  }
  #### Filter phenotypes ####
  if(!is.null(keep_tiers)){
    keep_ids <- HPOExplorer::hpo_tiers[tier %in% keep_tiers,]$hpo_id
    results <- results[HPO_ID %in% keep_ids,]
  }
  #### Filter celltypes ####
  if(!is.null(keep_celltypes)){
    results <- results[tolower(CellType) %in% tolower(keep_celltypes),]
  }
  #### Filter genes ####
  p2g <- phenotype_to_genes[
    tolower(Phenotype) %in% tolower(unique(results$Phenotype)),
  ]
  spec <- ctd[[annotLevel]]$specificity
  shared_genes <- intersect(p2g$Gene,rownames(spec))
  spec <- spec[shared_genes,]
  p2g <- p2g[Gene %in% shared_genes]

}

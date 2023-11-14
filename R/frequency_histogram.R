#' Gene and phenotype frequencies histogram
#'
#' Plot the frequency of gene-phenotype and phenotype-disease associations.
#' @inheritParams prioritise_targets
#' @inheritParams frequency_histogram
#' @inheritParams prioritise_targets_network
#' @return ggplot object
#'
#' @export
#' @importFrom HPOExplorer load_phenotype_to_genes
#' @importFrom HPOExplorer add_gene_frequency add_ancestor
#' @examples
#' results <- load_example_results()[seq(5000),]
#' fp_res <- frequency_histogram(results=results)
frequency_histogram <- function(results = load_example_results(),
                                phenotype_to_genes = load_phenotype_to_genes(),
                                show_plot = FALSE,
                                verbose = TRUE){

  Metric_Type <- value <- variable <- hpo_id <- NULL;

  results <- results |>
    HPOExplorer::add_pheno_frequency() |>
    HPOExplorer::add_ancestor()
  gene_df <- phenotype_to_genes[hpo_id %in% unique(results$hpo_id),] |>
    HPOExplorer::add_gene_frequency() |>
    HPOExplorer::add_ancestor()

  #### Gene freqs ####
  measure.vars <- grep(paste("^p$","^q$","^fold_change",
                             "_min$","_max$","_mean$", sep = "|"),
                       names(gene_df), value = TRUE)
  d1 <- data.table::melt.data.table(
    gene_df,
    id.vars = c("hpo_id","hpo_name"),
    measure.vars = measure.vars )

  g1 <- ggplot2::ggplot(d1, ggplot2::aes(x=value, fill=variable)) +
    ggplot2::geom_histogram(stat = "count", na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = pals::viridis(4)) +
    ggplot2::facet_wrap(facets = "variable ~.") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background =
                     ggplot2::element_rect(fill = "transparent"))

  #### Phenotype freqs ####
  measure.vars <- grep(paste("^p$","^q$","^fold_change",
                             "_min$","_max$","_mean$", sep = "|"),
                       names(results), value = TRUE)
  d2 <- data.table::melt.data.table(results,
                                         id.vars = c("hpo_id","hpo_name"),
                                         measure.vars = measure.vars)
  d2$Metric_Type <- ifelse(grepl("pheno_freq",d2$variable),
                                "Frequency","Enrichment")
  g2 <- ggplot2::ggplot(d2,
                  ggplot2::aes(x=value, fill=Metric_Type) )+
    ggplot2::geom_density(alpha=.75, color = NA) +
    ggplot2::scale_fill_manual(values = pals::inferno(3)) +
    ggplot2::facet_wrap(facets = "variable ~.",
                        scales = "free", ncol = 4) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background =
                     ggplot2::element_rect(fill = "transparent"))
  #### Merge plots ####
  pw <- patchwork::wrap_plots(list(g1, g2), ncol = 1, heights = c(.5,1)) +
    patchwork::plot_annotation(tag_levels = LETTERS)
  #### Show plot ####
  if(isTRUE(show_plot)) methods::show(pw)
  return(list(plot=pw,
              data=list(gene_plot=d1,
                        pheno_plot_df=d2))
         )
}

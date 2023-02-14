#' Gene and phenotype frequencies bar plot
#'
#' Plot the frequency of gene-phenotype and phenotype-disease associations.
#' @inheritParams prioritise_targets
#' @inheritParams prioritise_targets_network
#' @return ggplot object
#'
#' @export
#' @importFrom HPOExplorer load_phenotype_to_genes
#' @examples
#' fp_res <- frequency_barplot()
frequency_barplot <- function(results = load_example_results(),
                              phenotype_to_genes = load_phenotype_to_genes(),
                              show_plot = TRUE,
                              verbose = TRUE){

  label <- percent <- Frequency <- ID <- NULL;

  results <- results |>
    HPOExplorer::add_pheno_frequency() |>
    HPOExplorer::add_ancestor()
  gene_df <- phenotype_to_genes[ID %in% unique(results$HPO_ID),] |>
    HPOExplorer::add_gene_frequency() |>
    HPOExplorer::add_ancestor()

  #### barplot #####
  #### Gene frequencies ####
  gene_plt_df <- frequency_plot_prepare(df = gene_df,
                                        col="gene_freq_name")
  gene_plt <- frequency_plot_barplot(
    plt_df = gene_plt_df,
    remove_x_text = TRUE,
    title="Gene frequency within phenotypes")
  #### Phenotype frequencies ####
  pheno_plt_df <- frequency_plot_prepare(df = results,
                                         col="pheno_freq_mean")
  pheno_plt <- frequency_plot_barplot(
    plt_df = pheno_plt_df,
    direction = 1,
    title="Phenotype frequency within diseases")

  pw <- patchwork::wrap_plots(list("Gene frequencies"=gene_plt,
                             "Phenotype frequencies"=pheno_plt),
                        ncol = 1) +
    patchwork::plot_annotation(tag_levels = LETTERS)

  if(isTRUE(show_plot)) methods::show(pw)
  return(list(plot=pw,
              data=list(gene_plot=gene_plt_df,
                        pheno_plot_df=pheno_plt_df)))
}

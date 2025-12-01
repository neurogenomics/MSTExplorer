#' @title Plot functions
#'
#' @description
#' Functions to generate plots.
#' @param x_var Variable to plot on the x-axis.
#' @param y_var Variable to plot on the y-axis.
#' @param fill_var Variable to fill by.
#' @param facet_var Variable to facet by.
#' @param show_plot Print the plot to the console.
#' @param save_path Save the plot to a file.
#' Set to \code{NULL} to not save the plot.
#' @param y_lab y-axis label.
#' @param x_lab x-axis label.
#' @param title Title of the plot.
#' @param subtitle Subtitle of the plot.
#' @param ncol Number of facet columns.
#' @param nrow Number of facet rows for the plot.
#' @param height Height of the saved plot.
#' @param width Width of the saved plot.
#' @param heights Passed to \link[patchwork]{wrap_plots}.
#' @param subtitle_size Size of the plot subtitle.
#' @param phenotype_to_genes Phenotype to gene mapping from
#' \link[HPOExplorer]{load_phenotype_to_genes}.
#' @param dpi Resolution of image after rasterization.
#' @family plot_
#' @returns R object.
#' @name plot_
NULL

#' Plot graph of phenotypes in the RD EWCE results subsetted by Cell type
#'
#' This subsets the results, using other functions from this script,
#' by cell type,
#' q value, and fold change. It then plots a graph/ network of
#'  the phenotypes and
#' colors the nodes by their fold change or q value
#' (see the \code{color_var} param).
#'
#'
#' \emph{Developer notes:}
#' The ontologyPlot package did not have prebuilt options to create the heatmap,
#' so the colours are assigned manually to each phenotype in this function.
#' There must be a more efficient way to do this.
#' @param color_var Name of the variable (column) in \code{results}
#' that should be used to color the graph nodes.
#' @inheritParams ggnetwork_plot_full
#' @inheritParams HPOExplorer::make_phenos_dataframe
#' @inheritParams scales::col_numeric
#' @inheritParams ontologyPlot::onto_plot
#' @inheritDotParams ontologyPlot::onto_plot
#' @returns A \pkg{ontologyPlot} plot of the network of phenotypes
#'  in a subset of RD EWCE results.
#'
#' @export
#' @importFrom scales col_numeric
#' @importFrom ontologyPlot onto_plot
#' @examples
#' plt <- ontology_plot(cell_type="Amacrine cells")
ontology_plot <- function(cell_type,
                          results = load_example_results(),
                          color_var = c("fold_change","q","p"),
                          q_threshold = 0.0005,
                          fold_threshold = 1,
                          hpo = HPOExplorer::get_hpo(),
                          phenotype_to_genes =
                            HPOExplorer::load_phenotype_to_genes(),
                          palette = "Spectral",
                          shape = "rect",
                          ...
                          ){
  # templateR:::source_all()
  # templateR:::args2vars(ontology_plot)

  HPO_term_valid <- NULL;

  message("Generating ontology plot.")
  cells <- get_cell_ontology(cell_type = cell_type,
                             results = results,
                             q_threshold = q_threshold,
                             fold_threshold = fold_threshold,
                             phenotype_to_genes = phenotype_to_genes,
                             hpo = hpo)
  cells <- cells[HPO_term_valid==TRUE,]
  ### Check color_var ####
  color_var <- color_var[[1]]
  val_opts <- eval(formals(ontology_plot)$color_var)
  if(!color_var %in% val_opts){
    messager("WARNING: color_var must be one of:",
             paste("\n -",shQuote(val_opts),collapse = ""));
    color_func <- function(x) 0
  } else{
    color_func <- scales::col_numeric(palette = "Spectral",
                                      cells[[color_var]],
                                      reverse = color_var == "fold change")
  }
  #### Create plot ####
  plt <- ontologyPlot::onto_plot(ontology = hpo,
                                 terms = cells$HPO_ID,
                                 fillcolor = color_func(cells[[color_var]]),
                                 shape = shape,
                                 ...)
  return(plt)
}


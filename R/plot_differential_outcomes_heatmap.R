#' Plot differential outcomes: heatmap
#' @export
#' @inheritDotParams add_symptom_results
#' @examples
#' results <- load_example_results()
#' keep_descendants <- "Hypotonia" # HP:0001252
#' results2 <- HPOExplorer::filter_descendants(results,
#'                                             keep_descendants = keep_descendants)
#' results2 <- HPOExplorer::add_death(results2,
#'                                   allow.cartesian = TRUE,
#'                                   agg_by = c("disease_id","hpo_id"))
#' out <- plot_differential_outcomes_heatmap(results=results2)
plot_differential_outcomes_heatmap <- function(phenotypes = NULL,
                                               results = load_example_results(),
                                               celltype_col="cl_name",
                                               outcome_var="AgeOfDeath_score_mean",
                                               x_var="celltype_symptom",
                                               title="Differential outcomes by cell type",
                                               subtitle=paste0(
                                                 "Phenotype(s): ",
                                                 paste(phenotypes,
                                                       collapse = "; ")),
                                               fill_limits=NULL,
                                               print_phenotypes=TRUE,
                                               show_plot=TRUE,
                                               save_path=NULL,
                                               height=NULL,
                                               width=NULL,
                                               ...){
  hpo_id <- NULL;
  results <- data.table::copy(results)
  if(!outcome_var %in% names(results)){
    stopper("outcome_var ",shQuote(outcome_var)," not in results.")
  }
  if(!is.null(phenotypes)){
    hpo_ids <- HPOExplorer::map_phenotypes(phenotypes, to = "id")
    results <- results[hpo_id %in% hpo_ids,]
  }
  # if(isTRUE(extend_phenotypes)){
  #   annot <- HPOExplorer::load_phenotype_to_genes(3)
  #   disease_ids <- annot[hpo_id %in% hpo_ids]$disease_id |>unique()
  #   hpo_ids <- annot[disease_id %in% disease_ids]$hpo_id |>unique()
  # }
  results <- add_symptom_results(results = results,
                                 celltype_col=celltype_col,
                                 ...)
  dat <- results[!is.na(get(outcome_var)) &
                 !is.na(get(x_var))]|>
    data.table::setorderv(outcome_var)
  #### sort by mean AgeOfDeath_score_mean and make an ordered factor ####
  dat$disease_name <- factor(dat$disease_name,
                             levels = unique(dat$disease_name),
                             ordered = TRUE)
  dat_mean <- data.table::copy(dat)[,(outcome_var):=mean(get(outcome_var)),
                                    by=c(x_var), ]|>
    data.table::setorderv(outcome_var)
  dat_mean[[x_var]] <- factor(dat_mean[[x_var]],
                             levels = unique(dat_mean[[x_var]]),
                             ordered = TRUE)
  dat[[x_var]] <- factor(dat[[x_var]],
                             levels = unique(dat_mean[[x_var]]),
                             ordered = TRUE)
  #### Determine legend fill limits ####
  if(is.null(fill_limits)){
    fill_limits <- c(min(dat[[outcome_var]]),
                     max(dat[[outcome_var]]))
  }

  #### Create subtitle ####
  if(is.null(phenotypes)){
    if(isTRUE(print_phenotypes)){
      dat <- HPOExplorer::add_hpo_name(dat)
      subtitle <- paste0("Phenotype(s): ",
                         paste(unique(dat$hpo_name),
                               collapse = "; "))
    } else {
      subtitle <- paste0("Phenotype(s): ",data.table::uniqueN(dat$hpo_id))
    }
  }
  #### Cell types means barplot ####
  g1 <- ggplot2::ggplot(dat_mean,
                        ggplot2::aes(x=!!ggplot2::sym(x_var),
                                     y=!!ggplot2::sym(outcome_var),
                                     fill=!!ggplot2::sym(outcome_var))) +
    ggplot2::geom_col(data = dat_mean[,.SD[1],by=c(x_var)],
                      show.legend = FALSE) +
    ggplot2::geom_boxplot(show.legend = FALSE)+
    ggplot2::scale_fill_viridis_c(option = "inferno", direction = -1,
                                  limits=fill_limits) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.margin = ggplot2::margin(0,0,0,0),
                   plot.background = ggplot2::element_blank()) +
    ### Shift over the ggplot grid lines so they are on either side of the bars
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::geom_vline(xintercept = seq(1,length(unique(dat_mean[[x_var]])))+.5,
                        alpha=.25, linewidth = 0.5)
  #### Disease x cell type heatmap ####
  g2 <- ggplot2::ggplot(dat,
                  ggplot2::aes(x=!!ggplot2::sym(x_var),
                               y=disease_name,
                               fill=!!ggplot2::sym(outcome_var))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "inferno", direction = -1,
                                  limits=fill_limits) +
    ggplot2::theme_bw() +
    ggplot2::labs(x=NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0,0,0,0),
                   plot.background = ggplot2::element_blank(),
                   legend.position = "top",
                   legend.direction = "horizontal") +
    ### Shift over the ggplot grid lines so they are on either side of the bars
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::geom_vline(xintercept = seq(1,length(unique(dat[[x_var]])))+.5,
                        alpha=.25, linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = seq(1,length(unique(dat[["disease_name"]])))+.5,
                        alpha=.25, linewidth = 0.5)
  #### Combine plots ####
  pw <- patchwork::wrap_plots(g2,
                        patchwork::plot_spacer(),
                        g1,
                        ncol=1,
                        # guides = "collect",
                        heights = c(1,-.1,.3)) +
    patchwork::plot_annotation(tag_levels = letters,
                               title=title,
                               subtitle=subtitle)
  #### Show ####
  if(isTRUE(show_plot)) methods::show(pw)
  #### Save ####
  KGExplorer::plot_save(plt  = pw,
                        path = save_path,
                        height = height,
                        width = width)
  #### Return ####
  return(
    list(plot=pw,
         data=dat)
  )
  # X <- data.table::dcast.data.table(long_data,
  #                                   formula=disease_name~cl_name,
  #                                   value.var = "AgeOfDeath_score_mean",
  #                                   fun.aggregate = mean,
  #                                   # fill = 0,
  #                                   na.rm=TRUE) |>
  #   as.data.frame()|>
  #   tibble::column_to_rownames("disease_name") |>
  #   as.matrix()
  # X <- X[names(sort(Matrix::rowMeans(X,na.rm=TRUE))),]
  # ComplexHeatmap::Heatmap(as.matrix(X),
  #                         cluster_columns = FALSE,
  #                         cluster_rows = FALSE)
}

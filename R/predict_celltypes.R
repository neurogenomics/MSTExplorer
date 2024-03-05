#' Predict cell types
#'
#' Predict the causal cell types underlying a patient's phenotypes given some
#' varying degree of prior knowledge.
#' @param phenotypes Phenotypes observed in the patient.
#' Can be a list of HPO phenotype IDs or HPO phenotype names.
#' @param diseases Diseases that the patient is known to have.
#' Can be provided as OMIM, Orphanet, or DECIPHER disease IDs.
#' @param genes_include Genes in which the patient is known to have
#'  abnormalities.
#' @param genes_exclude Genes in which the patient is known NOT to have
#'  abnormalities.
#' @inheritParams prioritise_targets
#' @returns \link[data.table]{data.table} of prioritised cell types,
#' sorted by a "score" that combines:
#' \itemize{
#' \item{The phenotype-cell type enrichment p-values ("p").}
#' \item{The phenotype-cell type enrichment effect size ("fold_change").}
#' \item{A gene-wise factor that upweights/downweights
#' included/excluded genes respectively,
#' multiplied by the evidence score of a phenotype-gene association.
#' Only applied when \code{genes_include} or \code{genes_exclude} is provided. }
#' }
#'
#' @export
#' @import data.table
#' @examples
#' phenotypes <- c("Generalized neonatal hypotonia",
#'                 "Scrotal hypospadias",
#'                 "Increased circulating progesterone")
#' # diseases_include <- "OMIM:176270"
#' genes_include <- c("MAGEL2","HERC2")
#' genes_exclude <- c("SNORD115-1")
#' ct <- predict_celltypes(phenotypes = phenotypes,
#'                         genes_include = genes_include,
#'                         genes_exclude = genes_exclude)
predict_celltypes <- function(phenotypes,
                              diseases_include=NULL,
                              diseases_exclude=NULL,
                              genes_include = NULL,
                              genes_exclude = NULL,
                              gene_weights = list(
                                include=2,
                                default=1,
                                exclude=0
                              ),
                              results = MSTExplorer::load_example_results(),
                              agg_var=c("cl_name"),#,"hpo_id","ctd"),
                              x_var=agg_var[1],
                              y_var="score_sum",
                              fill_var="score_mean",
                              max_x_var=10,
                              subtitle_size=9,
                              plot.margin = ggplot2::margin(1,1,1,40),
                              show_plot = TRUE,
                              save_path = NULL,
                              width=NULL,
                              height=NULL){
  requireNamespace("ggplot2")
  gene_symbol <- evidence_score_mean <- fold_change <- hpo_id <-
    disease_id <- p <- NULL;

  {
    hpo_ids <- HPOExplorer::map_phenotypes(terms = phenotypes,
                                           to = "id")
    res <- data.table::copy(results)
    # res <- res[q<0.05]
    # res[,fold_change:=scales::rescale(abs(fold_change))]
    #### Filter phenotypes ####
    if(!is.null(hpo_ids)) {
      res <- res[hpo_id %in% hpo_ids,]
    }
    #### Add diseases ####
    res <- HPOExplorer::add_disease(phenos = res,
                                    allow.cartesian = TRUE)
    #### Filter diseases ####
    if(!is.null(diseases_include)){
      res <- res[disease_id %in% diseases_include,]
    }
    if(!is.null(diseases_exclude)){
      res <- res[!disease_id %in% diseases_exclude,]
    }
    #### Add cell types ####
    res <- MSTExplorer::map_celltype(results = res)
  }
  #### Add genes ####
  if(!is.null(genes_include) ||
     !is.null(genes_exclude)){
    old_cols <- names(res)
    res <- HPOExplorer::add_genes(phenos = res)
    # res <- MSTExplorer:::add_driver_genes(res,
    #                                       # min_value = 0,
    #                                       metric = "specificity")
    # res <- MSTExplorer::add_symptom_results(res,
    #                                         q_threshold = 1,
    #                                         fold_threshold = 2,
    #                                         celltype_col = x_var,
    #                                         keep_quantiles = seq(30,40),
    #                                         proportion_driver_genes_symptom_threshold = .75)
    res <- HPOExplorer::add_evidence(phenos = res)
    res[,gene_factor:=ifelse(gene_symbol %in% genes_include,
                             gene_weights$include,
      ifelse(gene_symbol %in% genes_exclude,
             gene_weights$exclude,
             gene_weights$default)
      ) * evidence_score_mean]
    #### Compute gene-adjusted combined score ####
    res[,score:=scales::rescale((1-p)*fold_change*gene_factor,c(0,1))]
    if("intersection" %in% names(res)){
      res <- res[,-c("intersection")]
    }
    res <- res[,unique(c(old_cols,"score")), with=FALSE]|>
      unique()
  } else {
   #### Compute combined score ####
   res[,score:=(1-p)*fold_change]
  }
  #### Aggregate scores ####
  {
    agg_var <- agg_var[agg_var %in% names(res)]
    res_agg <- res[,list(p_mean=mean(p),
                         q_mean=mean(q),
                         fold_change_mean=mean(fold_change),
                         score_sum=sum(score),
                         score_mean=mean(score),
                         score_median=median(score),
                         score_sd=stats::sd(score)),
                   by=agg_var] |>
      data.table::setorderv(cols = "score_sum",
                            order = -1)
    res_agg[[x_var]] <- factor(res_agg[[x_var]],
                               unique(res_agg[[x_var]]),
                               ordered = TRUE)
    res[[x_var]] <- factor(res[[x_var]],
                           as.character(unique(res_agg[[x_var]])),
                           ordered = TRUE)
  }
  #### Create title (only including non-null arguments) ####
  {
    sep <- "; "
    title <- paste0("Phenotypes: ",
                    paste(phenotypes,collapse=sep))
    if(!is.null(diseases_include)){
      title <- paste0(title,"\nDiseases included: ",
                      paste(diseases_include,collapse=sep))
    }
    if(!is.null(diseases_exclude)){
      title <- paste0(title,"\nDiseases excluded: ",
                      paste(diseases_exclude,collapse=sep))
    }
    if(!is.null(genes_include)){
      title <- paste0(title,"\nGenes included: ",
                      paste(genes_include,collapse=sep))
    }
    if(!is.null(genes_exclude)){
      title <- paste0(title,"\nGenes excluded: ",
                      paste(genes_exclude,collapse=sep))
    }
  }
  #### Plot ####
  res[,score_mean:=mean(score), by=x_var]
  res[,score_sum:=sum(score), by=x_var]
  res[,score_median:=median(score), by=x_var]
  res[,n_phenotypes:=length(unique(hpo_id)), by=x_var]
  plot_dat <- if(!is.null(max_x_var)){
    res[get(x_var) %in% utils::head(levels(get(x_var)),max_x_var)]
  } else{
    data.table::copy(res)
  }
  gp <- ggplot2::ggplot(plot_dat,
                        ggplot2::aes(x=!!ggplot2::sym(x_var),
                                     y=!!ggplot2::sym(y_var),
                                     fill=!!ggplot2::sym(fill_var))) +
    ggplot2::geom_violin(show.legend = FALSE) +
    ggplot2::geom_col(data=plot_dat[,.SD[1],by=c(x_var)],
                      ggplot2::aes(y=!!ggplot2::sym(y_var))) +
    ggplot2::stat_summary(fun = "median",
                          geom = "crossbar",
                          # mapping = aes(color=score_sum),
                          show.legend = FALSE) +
    ggplot2::geom_point(alpha=.5, show.legend = FALSE) +
    ## 1 standard deviation above the mean
    ggplot2::geom_text(inherit.aes = FALSE,
                       check_overlap=TRUE,
                       y = (mean(res[[y_var]])+sd(res[[y_var]]))+max(res[[y_var]])*.03,
                       x=10,
                       nudge_y = 1,
                       hjust=1, size=3, alpha=.75,
                       label=paste0("mean(",y_var,") + 1 SD")) +
    ggplot2::geom_hline(yintercept = mean(res[[y_var]])+sd(res[[y_var]]),
                        linetype="dotted", alpha=.75) +
    ## The mean
    ggplot2::geom_text(inherit.aes = FALSE,
                       check_overlap=TRUE,
                       y = mean(res[[y_var]])+max(res[[y_var]])*.03,
                       x=10,
                       hjust=1, size=3, alpha=.75,
                       label=paste0("mean(",y_var,")")) +
    ggplot2::geom_hline(yintercept = mean(res[[y_var]]),
                        linetype="dashed", alpha=.75) +
    ggplot2::labs(subtitle = title) +
    # annotate(x = .75, y=mean(res$score)*1.05,
    #          geom = "text",
    #          label="global mean",
    #          size=3,
    #          alpha=.5) +
    ggplot2::scale_fill_gradient(low = "blue", high = "red") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust=1),
                   plot.margin = plot.margin,
                   plot.subtitle = ggplot2::element_text(size=subtitle_size))
  if(isTRUE(show_plot)) methods::show(gp)
  #### Save ####
  KGExplorer::plot_save(gp,
                        path = save_path,
                        width = width,
                        height = height)
  #### Return #####
  return(list(data=res,
              data_agg=res_agg,
              plot=gp))
}

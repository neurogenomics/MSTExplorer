#' Predict cell types
#'
#' Predict the causal cell types underlying a patient's phenotypes given some
#' varying degree of prior knowledge.
#' @param phenotypes Phenotypes observed in the patient.
#' Can be a list of HPO phenotype IDs or HPO phenotype names.
#' @param diseases_include Diseases that the patient is known to have.
#' Can be provided as OMIM, Orphanet, or DECIPHER disease IDs.
#' @param diseases_exclude Diseases that the patient is known NOT to have.
#' Can be provided as OMIM, Orphanet, or DECIPHER disease IDs.
#' @param genes_include Genes in which the patient is known to have
#'  abnormalities.
#' @param genes_exclude Genes in which the patient is known NOT to have
#'  abnormalities.
#' @param gene_weights A named list describing the weight to apply to
#' genes in the include, default, and exclude lists.
#' @param agg_var The variable(s) to aggregate \code{results} by.
#' @param max_x_var The maximum number of cell types to display.
#' @param evidence_score_var Which variable from
#' \link[HPOExplorer]{add_evidence} to use when weighting genes.
#' @inheritParams plot_
#' @inheritParams prioritise_targets
#' @inheritParams ggplot2::theme
#' @returns \link[data.table]{data.table} of prioritised cell types,
#' sorted by a "score" that combines:
#' \itemize{
#' \item{The phenotype-cell type enrichment p-values ("p").}
#' \item{The phenotype-cell type enrichment effect size ("effect").}
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
                              phenotype_to_genes = HPOExplorer::load_phenotype_to_genes(),
                              agg_var=c("cl_name"),#,"hpo_id","ctd"),
                              effect_var="logFC",
                              x_var=agg_var[1],
                              y_var="score_mean",
                              fill_var="score_sum",
                              evidence_score_var="evidence_score_sum",
                              max_x_var=10,
                              subtitle_size=9,
                              plot.margin = ggplot2::margin(1,1,1,40),
                              show_plot = TRUE,
                              save_path = NULL,
                              width=NULL,
                              height=NULL){
  requireNamespace("ggplot2")
  gene_symbol <- hpo_id <-
    disease_id <- p <- gene_factor <- score <- score_mean <-
    score_sum <- score_median <- n_phenotypes <- specificity <- NULL;

  {
    hpo_ids <- HPOExplorer::map_phenotypes(terms = phenotypes,
                                           to = "id")
    res <- data.table::copy(results)
    add_logfc(res)
    res[,effect:=scales::rescale(get(effect_var),c(0,1))]
    # res <- res[q<0.05]
    # res[,effect:=scales::rescale(abs(effect))]
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
    res <- add_driver_genes(res,
                            phenotype_to_genes=phenotype_to_genes,
                            metric = "specificity")
    res <- HPOExplorer::add_evidence(phenos = res)
    res[,gene_factor:=ifelse(gene_symbol %in% genes_include,
                             gene_weights$include,
      ifelse(gene_symbol %in% genes_exclude,
             gene_weights$exclude,
             gene_weights$default)
      ) * get(evidence_score_var)]
    # res[,gene_factor:=scales::rescale(gene_factor,c(0,1))]
    #### Compute gene-adjusted combined score ####
    res[,score:=scales::rescale((1-q)*effect*gene_factor*specificity,c(0,1))]
    if("intersection" %in% names(res)){
      res <- res[,-c("intersection")]
    }
    res <- res[,unique(c(old_cols,"score")), with=FALSE]|>
      unique()
  } else {
   #### Compute combined score ####
   res[,score:=(1-q)*effect]
  }
  #### Aggregate scores ####
  {
    agg_var <- agg_var[agg_var %in% names(res)]
    res_agg <- res[,list(p_mean=mean(p),
                         q_mean=mean(q),
                         effect_mean=mean(effect),
                         score_sum=sum(score),
                         score_mean=mean(score),
                         score_median=stats::median(score),
                         score_sd=stats::sd(score)),
                   by=agg_var] |>
      data.table::setorderv(cols = y_var,
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
  {
    res[,score_mean:=mean(score), by=x_var]
    res[,score_sum:=sum(score), by=x_var]
    res[,score_median:=stats::median(score), by=x_var]
    res[,n_phenotypes:=length(unique(hpo_id)), by=x_var]
    plot_dat <- if(!is.null(max_x_var)){
      res[get(x_var) %in% utils::head(levels(get(x_var)),max_x_var)]
    } else{
      data.table::copy(res)
    }
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
    ## 2 standard deviations above the mean
    ggplot2::geom_text(inherit.aes = FALSE,
                       check_overlap=TRUE,
                       y = (mean(res[[y_var]])+stats::sd(res[[y_var]])*2)+max(res[[y_var]])*.03,
                       x=10,
                       nudge_y = 1,
                       hjust=1, size=3, alpha=.75,
                       label=paste0("mean(",y_var,") + 2 SD")) +
    ggplot2::geom_hline(yintercept = mean(res[[y_var]])+stats::sd(res[[y_var]]*2),
                        linetype="dotted", alpha=.75) +
    ## 1 standard deviation above the mean
    ggplot2::geom_text(inherit.aes = FALSE,
                       check_overlap=TRUE,
                       y = (mean(res[[y_var]])+stats::sd(res[[y_var]]))+max(res[[y_var]])*.03,
                       x=10,
                       nudge_y = 1,
                       hjust=1, size=3, alpha=.75,
                       label=paste0("mean(",y_var,") + 1 SD")) +
    ggplot2::geom_hline(yintercept = mean(res[[y_var]])+stats::sd(res[[y_var]]),
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
    # ggside::geom_ysidedensity(data=res,  alpha=0.5)
  if(isTRUE(show_plot)) methods::show(gp)
  #### Save ####
  KGExplorer::plot_save(gp,
                        save_path = save_path,
                        width = width,
                        height = height)
  #### Return #####
  return(list(data=res,
              data_agg=res_agg,
              data_plot=plot_dat,
              plot=gp))
}

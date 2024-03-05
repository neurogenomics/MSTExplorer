#' Plot bar plot of enriched phenotypes per celltype per HPO branch
#'
#' @export
#' @examples
#' results <- load_example_results()
#' results <- HPOExplorer::add_ancestor(results,
#'                                      lvl = 7,
#'                                      force_new = TRUE)
#' target_branches <- list("Recurrent bacterial infections"="leukocyte")
#' out <- plot_bar_branches(results=results,
#'                          target_branches=target_branches,
#'                          facets = "hpo_name",
#'                          legend.position="right",
#'                          lvl=8,
#'                          ncol=2,
#'                          vbars="hepatoblast",
#'                          facets_n=NULL,
#'                          q_threshold=0.05,
#'                          background_full=FALSE)
plot_bar_branches <- function(results = load_example_results(),
                              results_full = NULL,
                              target_branches = get_target_branches(),
                              keep_ancestors = names(target_branches),
                              target_celltypes = get_target_celltypes(
                                target_branches=target_branches
                              ),
                              target_branches_keep = NULL,
                              celltype_col = "cl_name",
                              keep_ont_levels=NULL,
                              add_test_target_celltypes=TRUE,
                              color_map=NULL,
                              color_vector=NULL,
                              preferred_palettes = "gnuplot",
                              legend.position="none",
                              fill_var="ancestor_name",
                              facets = fill_var,
                              scales="free_y",
                              q_threshold=NULL,
                              lvl=NULL,
                              facets_n="phenotypes_per_ancestor",
                              suffix="phenotypes",
                              vbars=NULL,
                              ncol=1,
                              cols=NULL,
                              y_lab="Significant phenotypes",
                              normalise_by=NULL,
                              background_full=TRUE,
                              save_path=NULL,
                              width=NULL,
                              height=NULL){
  requireNamespace("ggplot2")
  hpo_id <- p.adj.signif <- ancestor_name_original <- ancestor_name <- NULL;

  results <- map_celltype(results)
  if(!"sig_phenotypes" %in% names(results)){
    results[, sig_phenotypes:=data.table::uniqueN(hpo_id[q<q_threshold],
                                                       na.rm = TRUE),
            by=c(celltype_col,"cl_id","ancestor","ancestor_name")]
  }
  if(!"phenotypes_per_ancestor" %in% names(results)){
    results[, phenotypes_per_ancestor:=data.table::uniqueN(hpo_id),
            by=c("ancestor","ancestor_name")]
  }
  if(!is.null(lvl)){
    if("ancestor_name" %in% names(results)){
      results[,ancestor_name_original:=ancestor_name]
    }
    results <- HPOExplorer::add_ancestor(results,
                                         lvl = lvl,
                                         force_new = TRUE)
    results[, phenotypes_per_ancestor:=data.table::uniqueN(hpo_id),
            by=c("ancestor","ancestor_name")]
  }
  if("hpo_name" %in% c(fill_var,facets)){
    results <- HPOExplorer::add_hpo_name(results)
  }
  #### Create filtered dataset for plotting ####
  {
    if(!is.null(keep_ancestors) &&
       "hpo_name" %in% names(results)){
      dat <- HPOExplorer::filter_descendants(phenos = results,
                                             keep_descendants = keep_ancestors)
    } else {
      dat <- data.table::copy(results)
    }
    if(nrow(dat)==0) stopper("0 associations remaining.")
    if(!is.null(keep_ont_levels) || "ontLvl" %in% c(facets,cols)){
      dat <- HPOExplorer::add_ont_lvl(dat,
                                      keep_ont_levels = keep_ont_levels)
    }
    if(!is.null(target_branches_keep)){
      dat <- dat[!get(facets) %in% names(target_branches)]
    }
  }

  if(!is.null(lvl)){
    #### Reassign target_celltypes by inheriting from ancestor_name_original ####
    new_ancestors <- unique(dat$ancestor_name)
    target_celltypes <- lapply(stats::setNames(new_ancestors,
                                               new_ancestors),
                               function(b){
      target_celltypes[
        unique(dat[ancestor_name==b,]$ancestor_name_original)
      ]|>unlist(use.names = FALSE) |> unique()
    })
  }
  #### Select background to test for enrichment agains t####
  results <- if(isTRUE(background_full)){
    if(is.null(results_full)){
      results_full <- data.table::copy(results)
    }
    results_full
  } else {
    dat
  }
  #### Run tests ####
  if(isTRUE(add_test_target_celltypes)){
    target_tests <- test_target_celltypes(results=results,
                                          tests="across_branches_per_celltype",
                                          target_celltypes = target_celltypes,
                                          q_threshold = q_threshold)
    if(nrow(target_tests[[1]])==0) {
      messager("0 tests returned. Skipping annotation.")
      add_test_target_celltypes <- FALSE
    } else {
      dat <- merge(dat,
                   target_tests[[1]][,-c("p")],
                   all.x = TRUE,
                   by=c("ancestor_name","cl_id"))
      dat[p.adj.signif=="ns",p.adj.signif:=NA]
    }
  }
  #### Filter data ####
  if(!is.null(q_threshold)){
    dat <- dat[q<q_threshold]
  }
  #### Make facets ordered ####
  dat[[facets]] <- factor(dat[[facets]],
                          # set facet var in the order of the fill var to avoid reordering
                          levels = rev(unique(names(target_celltypes))),
                          # levels=unique(dat[order(match(get(facets),get(fill_var))),][[facets]]),
                          ordered = TRUE)
  if(!is.null(normalise_by) && normalise_by %in% names(dat)){
    dat[,sig_phenotypes:=scales::rescale_max(sig_phenotypes),
        by=normalise_by]
  }
  #### Create color map ####
  if(is.null(color_map)){
    color_map <- KGExplorer::map_colors(dat,
                                        columns = fill_var,
                                        preferred_palettes = preferred_palettes,
                                        as="dict")[[1]]
  }
  #### Bar plot ####
  ggbars <- ggplot2::ggplot(dat,
                            ggplot2::aes(x=!!ggplot2::sym(celltype_col),
                                         y=sig_phenotypes,
                                         fill=!!ggplot2::sym(fill_var)
                                         )
                            ) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::labs(x=NULL,
                  y=y_lab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = legend.position,
      strip.background = ggplot2::element_rect(fill = "transparent"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90,
                                          hjust = 1,
                                          vjust = 0.5)
    )

  #### Add facets ####
  if(!is.null(cols)){
    ggbars <- ggbars +
      ggplot2::facet_grid(facets = paste(facets,"~",cols),
                          scales = scales,
                          labeller = construct_labeller(dat=dat,
                                                        facets=facets,
                                                        facets_n=facets_n,
                                                        suffix=suffix))
  } else {
    ggbars <- ggbars +
      ggplot2::facet_wrap(facets =facets,
                          as.table = FALSE,
                          scales = scales,
                          labeller = construct_labeller(dat=dat,
                                                        facets=facets,
                                                        facets_n=facets_n,
                                                        suffix=suffix),
                          ncol = ncol)
  }

  if(!is.null(vbars)){
    ggbars <- ggbars +
      ggplot2::geom_vline(xintercept = vbars,
                          color="grey",alpha=1, linetype="dashed")
  }
  if(!is.null(color_map)){
    ggbars <- ggbars +
      ggplot2::scale_fill_manual(values=color_map)
  }
  if(!is.null(color_vector)){
    ggbars <- suppressWarnings(
      ggbars + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust=0.5,
                                            color = unname(color_vector))
      )
    )
  }
  #### Add test results to bar plot ####
  if(isTRUE(add_test_target_celltypes)){
    ggbars <- ggbars +
      ggplot2::geom_text(ggplot2::aes(label=p.adj.signif,
                                      y=1.05*sig_phenotypes),
                                      # nudge_y=3,
                                      size=2,
                                      color="black",
                                      alpha=.8,
                                      na.rm = TRUE) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .2)))
  }
  #### Save plot ####
  if(!is.null(save_path)){
    KGExplorer::plot_save(plt = ggbars,
                          path = save_path,
                          height = height,
                          width = width)
  }
  return(list(plot=ggbars,
              dat=dat))
}

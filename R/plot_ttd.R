#' Plot TTD
#'
#' Plot overlap between your targets and TTD targets.
#' @inheritParams ggplot2::labs
#' @keywords internal
plot_ttd <- function(dat_sub,
                     base_size = 9,
                     label_size = 2,
                     as_percent=TRUE,
                     baseline=NULL){
  requireNamespace("ggplot2")
  requireNamespace("scales")
  HIGHEST_STATUS<- prioritised <- n_drugs <- failed <- baseline_prop <- NULL;

  dat_sub[,n_drugs:=.N,
          by=c("HIGHEST_STATUS","prioritised")]
  dat_sub[,n_drugs_total:=.N, by=c("HIGHEST_STATUS")]
  ## Make function
  plot_ttd_make <- function(dat_sub,
                            subtitle_prefix=NULL,
                        subtitle = paste0(
                          subtitle_prefix,
                          "(",
                            round(sum(dat_sub[prioritised==TRUE,]$n_drugs)/
                                    sum(dat_sub$n_drugs)*100,2),
                            "% prioritised)"
                          ),
                        width = NULL,
                        show.legend=TRUE,
                        as_percent=TRUE,
                        ylims=NULL){

    if(isTRUE(as_percent)){
      plt <- ggplot2::ggplot(dat_sub,
                      ggplot2::aes(x=HIGHEST_STATUS,
                                   fill=prioritised,
                                   label=n_drugs)) +
        ggplot2::geom_bar(position = ggplot2::position_fill(),
          width = width,
          show.legend = show.legend) +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::geom_text(
          data = unique(dat_sub[,c("HIGHEST_STATUS","n_drugs_total")]),
          size=label_size,
          angle=90,
          check_overlap = TRUE,
          ggplot2::aes(x=HIGHEST_STATUS,
                       y=.95,
                       label=paste0("n=",n_drugs_total),
                       ),
          hjust=1,
          inherit.aes = FALSE)
      y_lab <- "Percent of therapeutic targets"
    } else {
      plt_dat <- dat_sub[,list(n_drugs),
                         by=c("HIGHEST_STATUS","prioritised","n_drugs_total"),
                         drop=FALSE] |> unique()
      plt_dat[prioritised==TRUE,
              prioritised_pct:=round((n_drugs/n_drugs_total*100),1),
              by=c("HIGHEST_STATUS")]

      plt <- ggplot2::ggplot(plt_dat,
                      ggplot2::aes(x=HIGHEST_STATUS,
                                   y=n_drugs,
                                   fill=prioritised,
                                   label=n_drugs)) +
        ggplot2::geom_col(width = width,
                          show.legend = show.legend) +
      ggplot2::geom_text(
        data = (data.table::setorderv(plt_dat,"prioritised_pct", na.last = TRUE)
                )[,.SD[1], by=c("HIGHEST_STATUS")],
        size=label_size,
        angle=90,
        check_overlap = TRUE,
        ggplot2::aes(label=paste0("  n=",n_drugs_total,
                                  "\n(",ifelse(is.na(prioritised_pct),0,prioritised_pct),
                                  "% prioritised)"),
                     x=HIGHEST_STATUS,
                     y=n_drugs_total
                     ),
        hjust= -.15,
        inherit.aes = FALSE) +
        ggplot2::scale_y_continuous(limits = ylims,
                                    expand = ggplot2::expansion(mult = c(0,.2))
                                    )
      ## Add baseline lines
      if(!is.null(baseline)){
        baseline_dt <- merge(baseline,
                             unique(plt_dat[,c("HIGHEST_STATUS","n_drugs_total")]),
                             by = "HIGHEST_STATUS")
        baseline_dt[,baseline_n:=(baseline_prop*n_drugs_total)]

        plt <- plt +
          ggplot2::geom_errorbar(data = baseline_dt,
                               ggplot2::aes(ymin = baseline_n,
                                            ymax = baseline_n,
                                            x = HIGHEST_STATUS),
                               linetype = 3,
                               linewidth = 1,
                               inherit.aes = FALSE,
                               show.legend = FALSE)
      }
      y_lab <- "Number of therapeutic targets"
    }
    plt <- plt +
      ggplot2::scale_fill_manual(guide = ggplot2::guide_legend(reverse=TRUE),
                                 values = c("grey80","springgreen3")) +
      ggplot2::labs(x="Approval stage",
                    y=y_lab,
                    fill="Prioritised",
                    subtitle=subtitle) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         hjust = 1))
    return(plt)
  }

  ## Make plot
  if("failed" %in% names(dat_sub)){
    fail <- dat_sub[failed==TRUE,drop=FALSE]
    notfail <- dat_sub[failed==FALSE,drop=FALSE]
    p1 <- plot_ttd_make(notfail,
                        subtitle_prefix="Non-failed therapeutics\n",
                        as_percent = as_percent)
    ylims <- c(0,
               max(p1$data$n_drugs_total))
    p2 <- plot_ttd_make(fail,
                        subtitle_prefix="Failed therapeutics\n",
                        as_percent = as_percent,
                        show.legend = FALSE,
                        ylims = ylims)
    prop_width <- length(unique(fail$HIGHEST_STATUS))/
                  length(unique(notfail$HIGHEST_STATUS))
    plt <- (p1 | p2) +
      patchwork::plot_layout(widths = c(1,prop_width),
                             axis_titles = "collect",
                             guides="collect",
                             axes="collect"
                             ) +
      patchwork::plot_annotation(tag_levels = "a")
  } else {
    plt <- plot_ttd_make(dat_sub,
                         as_percent = as_percent)
  }
  return(plt)
}

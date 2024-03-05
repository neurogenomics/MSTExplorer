#' Plot TTD
#'
#' Plot overlap between your targets and TTD targets.
#' @inheritParams ggplot2::labs
#' @keywords internal
plot_ttd <- function(dat_sub,
                     label_size=3,
                     failed_status=NULL){
  requireNamespace("ggplot2")
  requireNamespace("scales")
  HIGHEST_STATUS<- prioritised <- n_drugs <- NULL;

  dat_sub[,n_drugs:=.N, by=c("HIGHEST_STATUS")]
  plot_ttd_make <- function(dat_sub,
                        subtitle=paste0(
                          "Prioritised: ",
                          round(sum(dat_sub$prioritised)/
                                nrow(dat_sub)*100,2),"%"),
                        width = NULL,
                        show.legend=TRUE){
    ggplot2::ggplot(dat_sub,
                   ggplot2::aes(x=HIGHEST_STATUS,
                                fill=prioritised,
                                label=n_drugs)) +
      ggplot2::geom_bar(position = ggplot2::position_fill(),
                        width = width,
                        show.legend = show.legend) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::scale_fill_manual(guide = ggplot2::guide_legend(reverse=TRUE),
                                 values = c("grey","springgreen3")) +
      ggplot2::geom_text(
        data = unique(dat_sub[,c("HIGHEST_STATUS","prioritised","n_drugs")]),
        size=label_size,
        angle=90,
        hjust=1,
        ggplot2::aes(x=HIGHEST_STATUS, fill=prioritised,
                     label=paste0("n=",n_drugs), y=.95),
        # fill="white", alpha=.75,
        inherit.aes = FALSE) +
      ggplot2::labs(x="Approval stage",
                    y="Percent TTD targets",
                    subtitle=subtitle) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  if(!is.null(failed_status)){
    fail <- dat_sub[HIGHEST_STATUS %in% failed_status,drop=FALSE]
    notfail <- dat_sub[!HIGHEST_STATUS %in% failed_status,drop=FALSE]

    p1 <- plot_ttd_make(fail, show.legend = FALSE)
    p2 <- plot_ttd_make(notfail)
    prop_width <- length(unique(fail$HIGHEST_STATUS))/
                  length(unique(notfail$HIGHEST_STATUS))
    plt <- patchwork::wrap_plots(p2,p1,nrow = 1,
                          widths = c(1,prop_width),
                          guides="collect") +
      patchwork::plot_layout(axis_titles = "collect",
                             axes="collect")
  } else {
    plt <- plot_ttd_make(dat_sub)
  }
  return(plt)
}

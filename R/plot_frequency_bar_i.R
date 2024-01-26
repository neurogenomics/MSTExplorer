plot_frequency_bar_i <- function(plt_df,
                                   y_lab = "Proportion of results",
                                   remove_x_text = FALSE,
                                   direction = -1,
                                   option = "inferno",
                                   title=NULL){
  label <- Percent <- Frequency <- NULL;

  gp <- ggplot2::ggplot(plt_df,
                        ggplot2::aes(x=label,
                                     y=Percent,
                                     fill=Frequency)) +
    ggplot2::geom_col(position="fill") +
    ggplot2::scale_fill_viridis_d(option = option,
                                  na.value = "grey",
                                  direction = direction,
                                  guide = if(direction==1){
                                    ggplot2::guide_legend(reverse = TRUE)
                                  }) +
    # ggplot2::geom_text(data = totals,
    #   ggplot2::aes(label = total, x = ancestor_name, y=100),
    #    position = "fill", inherit.aes = FALSE, angle = 45, hjust=0
    # ) +
    ggplot2::labs(y=y_lab,
                  x=if(isFALSE(remove_x_text))"ancestor_name"else NULL,
                  title=title) +
    ggplot2::theme_bw()
  if(isFALSE(remove_x_text)){
    gp <- gp + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                                  hjust = 1))
  } else {
    gp <- gp + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }

  return(gp)
}

ttd_plot <- function(dat_sub){

  requireNamespace("ggplot2")
  HIGHEST_STATUS<- prioritised <- n_drugs <- NULL;

  dat_sub[,n_drugs:=.N, by=c("HIGHEST_STATUS")]
  ggplot(dat_sub, aes(x=HIGHEST_STATUS,
                      fill=prioritised,
                      label=n_drugs)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    geom_label(
      data = unique(dat_sub[,c("HIGHEST_STATUS","prioritised","n_drugs")]),
      aes(x=HIGHEST_STATUS, fill=prioritised, label=n_drugs, y=1),
      fill="white", alpha=.75, inherit.aes = FALSE) +
    labs(x="Approval stage") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

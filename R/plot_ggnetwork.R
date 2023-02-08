plot_ggnetwork <- function(g){

  g2 <- ggnetwork::ggnetwork(g)

  gp <- ggplot2::ggplot(g2,
                  ggplot2::aes(x=x, y=y,
                               xend = xend, yend = yend,
                               shape = node_type,
                               fill = node_type,
                               label = name)) +
    ggnetwork::geom_edges(ggplot2::aes(color = group),
                          curvature = 0.1,
                          alpha=0.25,
                          color = "grey50") +
    ggplot2::scale_fill_viridis_d(option = "magma",begin = .4) +
    # ggplot2::scale_shape_manual(values = shape_dict) +
    ggnetwork::geom_nodes(size=10) +
    # ggnetwork::geom_nodelabel(size=2, alpha=.8) +
    ggnetwork::geom_nodetext(size=2, alpha=.8) +
    ggnetwork::theme_blank()

  methods::show(gp)
  plotly::ggplotly(gp)
  # plotly::plot_ly(x=g2$x,
  #                 y=g2$y,
  #                 type = "scatter")
}

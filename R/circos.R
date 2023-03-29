## #' @export
## #' @import ggplot2
## #' @importFrom stats setNames
## #' @importFrom dplyr %>%
## #' @examples
## #' res <- MultiEWCE::example_targets
## #' vn <- prioritise_targets_network(top_targets = res$top_targets)
## prioritise_targets_circos <- function(top_targets){
##   # devoptera::args2vars(prioritise_targets_network, reassign = TRUE)
##   # http://blog.schochastics.net/post/visualizing-multilevel-networks-with-graphlayouts/
##   # https://mr.schochastics.net/material/netVizR/
##
##   g <- targets_to_graph(top_targets = top_targets,
##                         vertex_vars = c(group_var,vertex_vars),
##                         group_var = group_var,
##                         edge_color_var = edge_color_var,
##                         edge_size_var = edge_size_var,
##                         mediator_var = mediator_var,
##                         node_palette = node_palette,
##                         verbose = verbose)
##
##   igraph::V(g)$lvl <- as.numeric(igraph::V(g)$node_type)
##
##   visNetwork::visIgraph(igraph = g, layout = "dendrogram")
##   library(ggraph)
##   # ggraph(pa,layout="stress")+
##   #   geom_edge_link0(width=0.2,colour="grey")+
##   #   geom_node_point(col="black",size=0.3)+
##   #   theme_graph()
##   bc <- betweenness(g)
##
##   ggraph(g,layout = "centrality", centrality = bc, tseq = seq(0,1,0.15)) +
##     draw_circle(use = "cent") +
##     annotate_circle(bc,format="",pos="bottom") +
##     geom_edge_link0(edge_color="black",edge_width=0.3)+
##     geom_node_point(aes(fill=as.factor(Faction)),size=2,shape=21)+
##     scale_fill_manual(values=c("#8B2323", "#EEAD0E"))+
##     theme_graph()+
##     theme(legend.position = "none")+
##     coord_fixed()+
##     labs(title="betweenness centrality")
##
## }

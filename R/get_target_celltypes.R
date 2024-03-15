get_target_celltypes <- function(target_branches = get_target_branches(),
                                 cl = KGExplorer::get_ontology("cl",
                                                               remove_rings=TRUE)
                                 ){
  KGExplorer::get_ontology_descendants(ont = cl,
                                       terms = target_branches)
}

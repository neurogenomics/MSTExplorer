get_target_celltypes <- function(target_branches = get_target_branches(),
                                 cl = get_cl()
                                 ){
  KGExplorer::get_ontology_descendants(ont = cl,
                                       terms = target_branches)
}

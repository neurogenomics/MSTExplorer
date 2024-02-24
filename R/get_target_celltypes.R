get_target_celltypes <- function(target_branches=get_target_branches(),
                                 cl = KGExplorer::get_ontology("cl")){
  lapply(target_branches, function(x){
    xt <- KGExplorer::map_ontology_terms(ont = cl,
                                         terms = x,
                                         to = 'id')
    simona::dag_offspring(cl,
                          include_self = TRUE,
                          term=xt)
  })
}

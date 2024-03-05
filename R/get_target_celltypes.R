get_target_celltypes <- function(target_branches=get_target_branches(),
                                 cl = KGExplorer::get_ontology("cl",
                                                               remove_rings=TRUE)
                                 ){
  lapply(target_branches, function(x){
    message(x)
    xt <- KGExplorer::map_ontology_terms(ont = cl,
                                         terms = x,
                                         to = 'id')
    if(all(is.na(xt))) {
      messager("WARNING: The term",x,"was not found in the ontology.")
      return(NULL)
    }
    simona::dag_offspring(cl,
                          include_self = TRUE,
                          term=xt)
  })
}

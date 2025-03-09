fix_cl_ids <- function(dat,
                       ont=get_cl(),
                       replace_map=list(
                         "CSH1_CSH2_positive_cells"="CL:0000351",
                         "SLC26A4_PAEP_positive_cells"="CL:0002097",
                         "Metanephric_cells"="CL:0002525",
                         "AFP_ALB_positive_cells"="CL:0005026",
                         "Parietal_and_chief_cells"="CL:0002180",
                         "PDE1C_ACSM3_positive_cells"="CL:1000313"
                       ),
                       obsolete_map=list(
                         "CL:0000003"="CL:0000000",
                         "CL:0000775"="CL:0000096"
                       )
                       ){
  # dat <- celltype_maps
  #### Find obsolete terms ####
  # obsolete_terms <- as.character(unique(dat$cl_id[!dat$cl_id %in% ont@terms]))
  #### Find unnecessarily vague terms ####
  # vague_terms <- unique(
  #   dat[cl_id %in% c( "CL:0000003","CL:0000000","~~all~~")]$author_celltype
  # )
  cl_name <- cl_id <- author_celltype <- NULL;
  if(!"cl_id" %in% names(dat)) stopper("No 'cl_id' column in the input data.")
  #### Replace generic celltype annotations with more specific ones ####
  if(length(replace_map)>0){
    for(x in names(replace_map)){
      if(nrow(dat[author_celltype==x])>0){
        dat <- dat[author_celltype==x, cl_id:=replace_map[author_celltype]]
      }
    }
  }
  if(length(obsolete_map)>0 &&
     nrow(dat[cl_id %in% names(obsolete_map)])>0){
    for(x in names(obsolete_map)){
      if(nrow(dat[cl_id==x])>0){
        dat[cl_id==x, cl_id:=obsolete_map[x]]
      }
    }
  }
  ### Relabel the cl_name based on the updated cl_id col
  dat$cl_name <- KGExplorer::map_ontology_terms(ont = ont,
                                 terms = as.character(dat$cl_id),
                                 to = 'name')
  unmapped_terms <- unique(dat[is.na(cl_name)]$cl_id)
  if(length(unmapped_terms)>0) {
    messager(length(unmapped_terms),
             " cl_name were not mapped to the ontology.")
  }
  return(dat)
}

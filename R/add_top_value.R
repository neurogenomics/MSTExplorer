add_top_value <- function(dat,
                          sort_var="sig_phenotypes",
                          label_var="ancestor_name",
                          group_var="cl_name",
                          new_var="top_ancestor_name",
                          normalise_group=TRUE
                          ){
  ## Normalise count within each ancestor
  if(normalise_group){
    sort_var_norm <- paste0(sort_var,"_norm")
    dat[,(sort_var_norm):=scales::rescale_max(get(sort_var)),
        by=label_var]
    sort_var <- sort_var_norm
  }
  dat[,(new_var):=head(get(label_var)[which.max(get(sort_var))],1),
      by=group_var]
}

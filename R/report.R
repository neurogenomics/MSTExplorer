report <- function(dt,
                   rep_dt = NULL,
                   step = NULL,
                   verbose = TRUE){

  if(!methods::is(dt,"data.table")){
    dt <- data.table::as.data.table(dt)
  }

  dict <- list(Rows="Rows",
               hpo_id="Phenotypes",
               disease_id="Diseases",
               CellType="Cell types",
               gene_symbol="Genes")
  get_values <- function(add_text=TRUE){
    lapply(stats::setNames(names(dict),names(dict)),
           function(nm){
             if(nm %in% names(dt) | (nm=="Rows")){
               if(nm=="Rows"){
                 val <- nrow(dt)
               } else {
                 val <- length(unique(dt[[nm]]))
               }
               if(add_text){
                 paste0(dict[[nm]],": ",formatC(val,big.mark = ","))
               } else{
                 val
               }
             } else {
               NULL
             }
           }) |> unlist()
  }
  rep <- get_values(add_text = TRUE)
  vals <- get_values(add_text = FALSE)
  messager("Prioritised targets:",
           if(!is.null(step))paste0("step=",shQuote(step)),
           paste("\n -",rep, collapse = " "), v=verbose)
  ##### Construct final table ####
  if(!is.null(step)){
    rep_dt2 <- data.table::data.table(
      step=step,
      t(vals |> `names<-`( dict[names(vals)])),
      ids=list(unique(dt$hpo_id))
    )
    return(data.table::rbindlist(list(rep_dt,rep_dt2), fill = TRUE))
  }
}

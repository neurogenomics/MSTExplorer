ttd_import <- function(save_dir = tools::R_user_dir(package = "MultiEWCE",
                                                    which = "cache"),
                       run_map_genes = TRUE){

  requireNamespace("orthogene")
  requireNamespace("readxl")
  requireNamespace("tidyr")
  TTDDRUID <- value <- GENENAME <- GENENAME2 <- GENENAME3 <- GENENAME_mapped <-
    TARGNAME <- HIGHEST_STATUS <- NULL;

  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  domain <- "https://db.idrblab.net/ttd/sites/default/files/ttd_database"
  files <- c("P1-02-TTD_drug_download.txt",
             "P1-01-TTD_target_download.txt",
             "P1-07-Drug-TargetMapping.xlsx",
             "P1-05-Drug_disease.txt")
  f <- lapply(stats::setNames(files,files),
              function(x){
    loc <- file.path(save_dir,x)
    if(!file.exists(loc)){
      utils::download.file(paste(domain,x,sep="/"),
                           destfile = loc)
    }
    loc
  })
  ##### Parse drug data #####
  drugs <- data.table::fread(f[["P1-02-TTD_drug_download.txt"]],
                             skip = 29,
                             col.names = c("DRUGID","key","value")) |>
    data.table::dcast.data.table(formula = "DRUGID  ~ key",
                                 value.var = "value",
                                 fun.aggregate = paste,
                                 collapse=";")
  #### Parse target data ####
  targets <- data.table::fread(f[["P1-01-TTD_target_download.txt"]],
                               skip=40,
                               col.names = c("TARGETID","key",
                                             "value","value2","status"))
  targets2 <- targets[key %in% c("GENENAME","TARGNAME","TARGTYPE","status"),] |>
    data.table::dcast.data.table(formula = "TARGETID  ~ key",
                                 value.var = "value",
                                 fun.aggregate = paste,
                                 collapse=";")
  #### Parse drug-target links ####
  drug_target <- readxl::read_excel(f[["P1-07-Drug-TargetMapping.xlsx"]])
  drug_target <- data.table::data.table(drug_target) |>
    `names<-`(toupper(names(drug_target)))
  #### Parse drug-disease data ####
  l <- suppressWarnings(
    readLines(f[["P1-05-Drug_disease.txt"]])
  )
  diseases <- data.table::fread(
    text = l[-seq(1,utils::tail(grep("----",l),1)+1)],
    col.names = c("key","value"))[key!="",]
  diseases <- rbind(list("TTDDRUID","D00ABE"),diseases)
  diseases[,TTDDRUID:=ifelse(key=="TTDDRUID",value,NA)]
  diseases <- tidyr::fill(diseases,TTDDRUID,.direction = "down")
  diseases2 <- data.table::dcast.data.table(diseases[key!="TTDDRUID",],
                                            formula = "TTDDRUID ~ key",
                                            value.var = "value",
                                            fun.aggregate = paste,
                                            collapse=";") |>
    data.table::setnames("TTDDRUID","DRUGID")
  ##### Merge #####
  dat <- drug_target[
    targets2, on="TARGETID"
  ][
    drugs, on="DRUGID"
  ][
    diseases2, on="DRUGID"
  ][,GENENAME2:=ifelse(GENENAME=="" | is.na(GENENAME),
                       stringr::str_extract(TARGNAME,
                                            pattern = "(?<=\\().*(?=\\))"),
                       GENENAME)][,GENENAME2:=strsplit(GENENAME2,"; ")] |>
    tidyr::unnest(cols = "GENENAME2") |>
    data.table::data.table()
  #### Harmonise gene names ####
  if(isTRUE(run_map_genes)){
    dat$GENENAME_mapped <- orthogene::map_genes(dat$GENENAME2,
                                                mthreshold = 1,
                                                drop_na = FALSE)$name
    dat[,GENENAME3:=dplyr::coalesce(GENENAME2,GENENAME_mapped)]
  }
  #### Convert status to ordered factor ####
  sts <- sort(unique(dat$HIGHEST_STATUS), na.last = TRUE)
  dat[,HIGHEST_STATUS:=stringr::str_to_sentence(HIGHEST_STATUS)]
  opts <- lapply(
    c(NA,
      "submitted|registered|preregistration|filed|patented",
      "investigative|preclinical",
      "phase|clinical",
      "approved",
      "discontinued|terminated|withdrawn"), FUN=function(x){
             grep(x,sts,value = TRUE,ignore.case = TRUE)
           }) |> unlist()|> unique()
  opts <- c(opts, setdiff(opts, sts))
  dat[,HIGHEST_STATUS:=factor(HIGHEST_STATUS, levels = opts, ordered = TRUE)]
  #### Return ####
  return(
    list(merged=dat,
         drugs=drugs,
         targets=targets,
         drug_target=drug_target,
         diseases=diseases)
  )
}

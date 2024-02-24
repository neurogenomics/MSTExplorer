#' Add tissues
#'
#' Annotate the tissues that each cell-type is found in.
#' @param celltypes A vector of cell-types to map onto tissues.
#' @param map A \link[data.table]{data.table} containing the columns
#'  "celltype" and "tissue".
#' If \code{NULL}, will use a built-in mapping file that only applies to the
#' CellTypeDataset DescartesHuman
#'  (see \link[MSTExplorer]{load_example_ctd} for details).
#' @param collapse If not \code{NULL},
#' collapse rows with >1 tissue into one string:
#' e.g. \code{c("Heart" "Lung")}  --> \code{"Heart;Lung"}
#' @returns A list of tissues with the same length
#'  as the input \code{celltypes}.
#' @source
#' \code{
#' #### Prepare mapping file #####
#' URL <- paste0("https://atlas.fredhutch.org/data/bbi/descartes/human_gtex/",
#'               "downloads/data_summarize_fetus_data/df_cell.RDS")
#' meta <- readRDS(url(URL))
#' map <- data.table::data.table(
#'   stringr::str_split(unique(meta$Organ_cell_lineage),"-",n = 2,
#'                      simplify = TRUE) |>
#'     `colnames<-`(c("tissue","celltype"))
#' )
#' f <- file.path("~/Desktop/ewce/rare_disease_celltyping/data/",
#'                "DescartesHuman_celltype_mapping.csv")
#'                data.table::fwrite(map,f)
#' piggyback::pb_upload(file = f, repo = "neurogenomics/MSTExplorer")
#' }
#'
#' @export
#' @examples
#' results <- load_example_results()
#' tissues <- map_tissues(celltypes = results$CellType)
map_tissues <- function(celltypes = NULL,
                        map = NULL,
                        collapse = NULL
                        ){
  tissue <- celltype <- NULL;

  if(is.null(map)){
    map <- get_data("DescartesHuman_celltype_mapping.csv")
    map[,celltype:=EWCE::fix_celltype_names(celltype,make_unique = FALSE)]
  }
  #### Count celltypes per tissue ####
  # map[,list(n=length(unique(celltype))), keyby="tissue"]
  map_agg <- map[,list(tissues=list(unique(tissue))), keyby="celltype"]
  ### Return early ####
  if(is.null(celltypes)) return(map)
  #### Subset ####
  celltypes <- EWCE::fix_celltype_names(celltypes = celltypes,
                                        make_unique = FALSE)
  tissues <- map_agg[celltypes,]$tissues
  #### Collapse lists into strings ####
  if(!is.null(collapse)){
    tissues <- sapply(tissues,paste, collapse=collapse)
  }
  #### Return ####
  return(tissues)
}

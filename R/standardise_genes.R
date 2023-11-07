#' Standardise genes
#'
#' Standardise gene symbols to HGNC.
#' @param dat data.table
#' @param gene_col character column name.
#' @param fill_na logical. Fill NAs with original gene symbol.
#' @param verbose logical. Print messages.
#' @inheritDotParams orthogene::map_genes
#' @returns data.table
#'
#' @export
#' @import data.table
#' @import orthogene
#' @examples
#' \dontrun{
#' dat <- data.table(gene_symbol = c("BRCA1","BRCA2","BRCA3"))
#' dat2 <- standardise_genes(dat)
#' }
standardise_genes <- function(dat,
                              gene_col = "gene_symbol",
                              fill_na = TRUE,
                              verbose = TRUE,
                              ...){

  input <- name <- equal <- gene_symbol_standard <- gene_symbol <- NULL;
  requireNamespace("orthogene", quietly = TRUE)
  gene_map <- orthogene::map_genes(genes = unique(dat[[gene_col]]),
                                   drop_na = TRUE,
                                   mthreshold = 1,
                                   # ...
                                   ) |>
    data.table::data.table()

  # data.table::uniqueN(gene_map$input)
  ## Prioritise exact matches
  gene_map[,equal:=input==name] |>
    data.table::setorderv("equal",-1)
  ## Get 1 row per gene
  gene_map <- gene_map[, .SD[1], by="input"] |>
    data.table::setnames("name", "gene_symbol_standard")
  ## Report any differences
  messager(length(unique(gene_map[input!=name,]$input)),
          "gene symbols differ from their standardised name.",v=verbose)
  dat <- data.table::merge.data.table(
    x = dat,
    y = gene_map[,c("input","gene_symbol_standard")],
    by.x = gene_col,
    by.y = "input",
    all.x = TRUE)
  ## Fill NAs
  if(isTRUE(fill_na)){
    dat[,gene_symbol_standard:=data.table::fcoalesce(
      gene_symbol_standard, gene_symbol)
      ]
  }
  return(dat)
}

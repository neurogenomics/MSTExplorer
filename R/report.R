report <- function(dt,
                   verbose = TRUE){
  messager(
    "Prioritised targets:",
    paste("\n -",formatC(nrow(dt),big.mark = ","),
          "results",
    "\n -",formatC(length(unique(dt$Phenotype)),big.mark = ","),
          "phenotypes",
    "\n -",formatC(length(unique(dt$CellType)),big.mark = ","),
          "cell types",
    "\n -",formatC(length(unique(dt$DiseaseName)),big.mark = ","),
          "associated diseases",
    "\n -",formatC(length(unique(dt$Gene)),big.mark = ","),
          "genes"
       ),
    v=verbose)
}

ttd_hypergeo <- function(fail,
                         notfail,
                         top_targets,
                         bg=c("p2g","fail","notfail"),
                         p2g=HPOExplorer::load_phenotype_to_genes()){
  # https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/

  ## Test for over-representation (enrichment)
  # phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)
  overlap <- data.table::uniqueN(notfail[prioritised==TRUE]$GENENAME3)
  group2 <- data.table::uniqueN(notfail$GENENAME3)
  total <- data.table::uniqueN(c(p2g$gene_symbol,
                                 fail$GENENAME3,
                                 notfail$GENENAME3))
  group1 <- data.table::uniqueN(top_targets$gene_symbol)
  nonfailed_enrichment <- stats::phyper(overlap-1,
                                        group2,
                                        total-group2,
                                        group1,
                                        lower.tail= FALSE)
  messager("Non-failed gene targets enrichment p-value:",nonfailed_enrichment)
  ## Test for under-representation (depletion)
  # phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE)
  overlap <- data.table::uniqueN(fail[prioritised==TRUE]$GENENAME3)
  group2 <- data.table::uniqueN(fail$GENENAME3)
  total <- data.table::uniqueN(c(p2g$gene_symbol,
                                 fail$GENENAME3,
                                 notfail$GENENAME3))
  group1 <- data.table::uniqueN(top_targets$gene_symbol)
  failed_depletion <- stats::phyper(overlap,
                                    group2,
                                    total-group2,
                                    group1,
                                    lower.tail= TRUE)
  messager("Failed gene targets depletion p-value:",failed_depletion)
  return(
    list(
      nonfailed_enrichment=nonfailed_enrichment,
      failed_depletion=failed_depletion
    )
  )
}

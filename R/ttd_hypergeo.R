ttd_hypergeo<- function(
    fail,
    notfail,
    top_targets,
    bg = c("p2g","fail","notfail"),
    p2g = HPOExplorer::load_phenotype_to_genes(),
    return_long = TRUE
) {

  # --------------------------------------------------------------------
  # 1. Define background universe
  # --------------------------------------------------------------------
  universe <- unique(c(p2g$gene_symbol,
                       fail$GENENAME3,
                       notfail$GENENAME3))

  top_vec <- unique(top_targets$gene_symbol)

  # Helper: build contingency table & stats -----------------------------
  compute_metrics <- function(target_set, status_label) {

    target_genes <- unique(target_set$GENENAME3)

    # counts for 2×2 table
    TP <- sum(top_vec %in% target_genes)
    FP <- sum(top_vec %in% setdiff(universe, target_genes))
    FN <- sum(setdiff(target_genes, top_vec) %in% universe)
    TN <- sum(setdiff(universe, c(top_vec, target_genes)) %in% universe)

    # Hypergeometric p-value (overrepresentation)
    overlap <- TP
    group2 <- length(target_genes)
    total  <- length(universe)
    group1 <- length(top_vec)

    if (status_label == "nonfailed") {
      p_hyper <- stats::phyper(overlap - 1,
                               group2,
                               total - group2,
                               group1,
                               lower.tail = FALSE)
    } else {
      # Depletion test (under-representation)
      p_hyper <- stats::phyper(overlap,
                               group2,
                               total - group2,
                               group1,
                               lower.tail = TRUE)
    }

    # Effect-size metrics -----------------------------------------------
    OR  <- (TP * TN) / max(FP * FN, .Machine$double.eps)
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    PPV <- TP / (TP + FP)
    NPV <- TN / (TN + FN)
    FDR <- 1 - PPV

    data.table::data.table(
      status = status_label,
      overlap = TP,
      targets_prioritised = group1,
      targets_in_status = group2,
      universe = total,
      p = p_hyper,
      OR = OR,
      sensitivity = sensitivity,
      specificity = specificity,
      PPV = PPV,
      NPV = NPV,
      FDR = FDR,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN
    )
  }

  # Compute tables ------------------------------------------------------
  res_nonfailed <- compute_metrics(notfail, "nonfailed")
  res_failed    <- compute_metrics(fail, "failed")

  if (return_long) {
    return(rbind(res_nonfailed, res_failed))
  } else {
    return(list(nonfailed = res_nonfailed,
                failed = res_failed))
  }
}

# ttd_hypergeo <- function(fail,
#                          notfail,
#                          top_targets,
#                          bg=c("p2g","fail","notfail"),
#                          p2g=HPOExplorer::load_phenotype_to_genes()){
#   # https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
#
#   ## Test for over-representation (enrichment)
#   # phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)
#   overlap <- data.table::uniqueN(notfail[prioritised==TRUE]$GENENAME3)
#   group2 <- data.table::uniqueN(notfail$GENENAME3)
#   total <- data.table::uniqueN(c(p2g$gene_symbol,
#                                  fail$GENENAME3,
#                                  notfail$GENENAME3))
#   group1 <- data.table::uniqueN(top_targets$gene_symbol)
#   nonfailed_enrichment <- stats::phyper(overlap-1,
#                                         group2,
#                                         total-group2,
#                                         group1,
#                                         lower.tail= FALSE)
#   messager("Non-failed gene targets enrichment p-value:",nonfailed_enrichment)
#   ## Test for under-representation (depletion)
#   # phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE)
#   overlap <- data.table::uniqueN(fail[prioritised==TRUE]$GENENAME3)
#   group2 <- data.table::uniqueN(fail$GENENAME3)
#   total <- data.table::uniqueN(c(p2g$gene_symbol,
#                                  fail$GENENAME3,
#                                  notfail$GENENAME3))
#   group1 <- data.table::uniqueN(top_targets$gene_symbol)
#   failed_depletion <- stats::phyper(overlap,
#                                     group2,
#                                     total-group2,
#                                     group1,
#                                     lower.tail= TRUE)
#   messager("Failed gene targets depletion p-value:",failed_depletion)
#   return(
#     list(
#       nonfailed_enrichment=nonfailed_enrichment,
#       failed_depletion=failed_depletion
#     )
#   )
# }

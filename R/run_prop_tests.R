#' Run proportional enrichment tests
#'
#' Run a series of proportional enrichment tests on the results of a
#' phenotype to cell type association test.
#' @keywords internal
#' @examples
#' results <- load_example_results()
#' results <- HPOExplorer::add_ancestor(results)
run_prop_tests <- function(results,
                           branch_col="ancestor_name",
                           celltype_col="CellType",
                           func=list(
                             rstatix::fisher_test,
                             rstatix::prop_test
                           )[[1]],
                           alternative = "greater",
                           cores=NULL,
                           ...){
  sig_phenos <- total_phenos <- NULL;
  sig_counts <- results[,list(sig_phenos=sum(q<.05),
                              nonsig_phenos=sum(q>=0.05),
                              total_phenos=.N),
                        by=c(branch_col,celltype_col)]
  BPPARAM <- KGExplorer::set_cores(workers = cores)
  test_res <- BiocParallel::bplapply(
    unique(sig_counts[[branch_col]]),
    BPPARAM = BPPARAM,
    function(anc){
     lapply(
      unique(sig_counts[[celltype_col]]),
      function(ct){
        res_null <- data.table::data.table()
        dt <- sig_counts[get(celltype_col)==ct][,target_ancestor:=get(branch_col)==anc]
        dt <- dt[,prop_phenos:=sig_phenos/total_phenos]
        if(all(dt$prop_phenos==0)) return (res_null)
        res <- tryCatch({
            xtab <- as.table(
              cbind(
                rbind(target_ancestor=sum(dt[target_ancestor==TRUE]$sig_phenos),
                      nontarget_ancestors=sum(dt[target_ancestor==FALSE]$sig_phenos)
                ),
                rbind(target_ancestor=sum(dt[target_ancestor==TRUE]$nonsig_phenos),
                      nontarget_ancestors=sum(dt[target_ancestor==FALSE]$nonsig_phenos)
                )
              )
            )
            dimnames(xtab)[[2]] <- c("sig_phenos","nonsig_phenos")
            func(xtab,
                 alternative = alternative,
                 ...)
        },
        error=function(e){print(e);return(res_null)})
        return(
          cbind(stats::setNames(list(ct),celltype_col),
                stats::setNames(list(anc),branch_col),
                N=dt[target_ancestor==TRUE]$sig_phenos,
                res)
        )
      })|> data.table::rbindlist(fill = TRUE)
    })|> data.table::rbindlist(fill = TRUE)
  # MTC
  if(nrow(test_res)>0){
    requireNamespace("gtools")
    test_res <- test_res |>
      dplyr::mutate(q=stats::p.adjust(p,"fdr"))|>
      dplyr::mutate(q_signif=gtools::stars.pval(p.value = q))
  }
  return(test_res)
}

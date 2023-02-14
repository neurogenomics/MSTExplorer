frequency_plot_prepare <- function(df,
                                   col="gene_freq_name"){
  ancestor_name <- n <- total <- pheno_freq <- NULL;
  requireNamespace("dplyr")

  plt_df <- data.table::copy(df)
  messager(paste0(
    formatC(round(sum(is.na(df[[col]]))/nrow(df)*100,1),big.mark = ","),"%"
  ),"of rows are NA.")
  if(grepl("^gene_freq",col)){
    get_freq_dict <- utils::getFromNamespace("get_freq_dict","HPOExplorer")
    plt_df[[col]] <- factor(plt_df[[col]],
                            levels = unname(get_freq_dict()),
                            ordered = TRUE)
  } else if(grepl("^pheno_freq",col)){
    plt_df <- plt_df[,pheno_freq:=mean(get(col), na.rm=TRUE),
                     by="ancestor_name"]
    plt_df$pheno_freq_og <- dplyr::ntile(plt_df$pheno_freq, n = 10)*10
    data.table::setorderv(plt_df,'pheno_freq_og',na.last = TRUE)
    lte <- "\U0002264"
    plt_df$pheno_freq <- gsub(paste0(lte,"NA%"),NA,
                              paste0(lte,plt_df$pheno_freq_og,"%"))
    plt_df$pheno_freq <- factor(plt_df$pheno_freq,
                               levels = unique(plt_df$pheno_freq),
                               ordered = TRUE)
    col <- "pheno_freq"
  }
  #### Count number of each frequency bin per ancestor_name ####
  plt_df <- plt_df |>
    dplyr::group_by_at(dplyr::all_of(c("ancestor_name",col))) |>
    dplyr::count() |>
    dplyr::rename(Frequency=eval(col))
  plt_df <- dplyr::group_by(plt_df, ancestor_name) |>
    dplyr::summarise(total=sum(n)) |>
    merge(x=plt_df, by="ancestor_name") |>
    #### Create x-axis labels with phenotype counts ####
    dplyr::mutate(Percent = n/total*100,
                  label=paste0(ancestor_name,"\n",
                               " (n=",formatC(total,big.mark = ","),")")
    )
  #### Order x-axis label #####
  plt_df$label <- factor(plt_df$label,
                         levels = unique(plt_df$label),
                         ordered = TRUE)

  return(plt_df)
}

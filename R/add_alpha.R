add_alpha <- function(col, alpha=.75){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, grDevices::col2rgb)/255, 2,
        function(x)
          grDevices::rgb(x[1], x[2], x[3], alpha=alpha))
}

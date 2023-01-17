#' Interactive DT
#'
#' Generate an interactive data table with download buttons.
#' @param dat Data to show.
#' @param caption Table caption.
#' @param scrollY Maximum height (in pixels) before defaulting to scrolling.
#'
#' @family general
#' @export
#' @examples
#' create_dt(dat = mtcars)
create_dt <- function(dat,
                      caption="",
                      scrollY=400){
  requireNamespace("DT")
    data <- DT::datatable(
      data = dat,
      caption = caption,
      extensions = 'Buttons',
      options = list( dom = 'Bfrtip',
                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                      scrollY = scrollY,
                      scrollX=TRUE, scrollCollapse = TRUE, paging = FALSE,
                      columnDefs = list(list(className = 'dt-center',
                                             targets = "_all"))
      )
    )
    return(data)
}

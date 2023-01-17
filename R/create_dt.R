#' Interactive DT
#'
#' Generate an interactive data table with download buttons.
#' @param dat Data to show.
#' @param caption Table caption.
#' @param scrollY Maximum height (in pixels) before defaulting to scrolling.
#' @param buttons Download button types to include in table.
#'
#' @family general
#' @export
#' @examples
#' create_dt(dat = mtcars)
create_dt <- function(dat,
                      caption = "",
                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                      scrollY = 400){
  requireNamespace("DT")
  DT::datatable(
    data = dat,
    caption = caption,
    extensions = 'Buttons',
    options = list( dom = 'Bfrtip',
                    buttons = buttons,
                    scrollY = scrollY,
                    scrollX=TRUE,
                    scrollCollapse = TRUE,
                    paging = FALSE,
                    columnDefs = list(list(className = 'dt-center',
                                           targets = "_all"))
    )
  )
}

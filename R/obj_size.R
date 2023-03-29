obj_size <- function(x,
                     units = "Mb",
                     digits = 5){
  format(utils::object.size(x),
         units = units,
         digits = digits)
}

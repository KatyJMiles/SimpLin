SimpLinR = function(x, y) {
  if (is.vector(x, mode = "numeric") && is.vector(y, mode = "numeric")) {
    SimpLinCpp(x,y)
  } else {
    stop("x and y need to be numeric vectors")
  }
}

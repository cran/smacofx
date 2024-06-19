#' Explicit Normalization
#' Normalizes distances
#' @param x numeric matrix 
#' @param w weight
#' @return a constant 
enorm <- function (x, w=1) {
    return (sqrt (sum (w * (x ^ 2))))
}

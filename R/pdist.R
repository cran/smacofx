#' Squared p-distances
#'
#' @param x numeric matrix
#' @param p p>0 the Minkoswki distance
#' @return squared Minkowski distance matrix
pdist <- function (x,p) {
    s <- tcrossprod (x)
    v <- diag (s)
    return (outer (v, v, "+") - 2 * s)
}

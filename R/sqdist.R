#' Squared distances
#'
#' @param x numeric matrix
#' @return squared distance matrix
sqdist <- function (x) {
    s <- tcrossprod (x)
    v <- diag (s)
    out <- outer (v, v, "+") - 2 * s #Was bug: could be less than 0
    out[out<0] <- 0 
    return(out)
}

sqdistOld <- function (x) {
    s <- tcrossprod (x)
    v <- diag (s)
    return(outer (v, v, "+") - 2 * s) #Bug
}

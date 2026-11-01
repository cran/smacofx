
#' Take matrix to a power 
#'
#' @param x matrix
#' @param r numeric (power)
#' @return a matrix
#' @noRd
mkPower<-function(x,r) {
    n<-nrow(x)
    tmp <- abs((x+diag(n))^r)-diag(n)
    tmp[!is.finite(tmp)] <- max(tmp[is.finite(tmp)]) #Was bug
    return(tmp)
}

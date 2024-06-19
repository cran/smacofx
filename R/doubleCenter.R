
#' Double centering of a matrix
#'
#' @param x numeric matrix
#' @return the double centered matrix
doubleCenter <- function(x) {
        n <- dim(x)[1]
        m <- dim(x)[2]
        s <- sum(x)/(n*m)
        xr <- rowSums(x)/m
        xc <- colSums(x)/n
        return((x-outer(xr,xc,"+"))+s)
    }

#'procruster: a procrustes function 
#'
#'@param x numeric matrix
#'@return a matrix
#'@noRd 
procruster <- function (x) 
{
    sx <- svd(x)
    return(tcrossprod(sx$u, sx$v))
}



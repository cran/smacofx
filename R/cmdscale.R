#' Wrapper to \code{cmdscale} for S3 class
#'
#' @details overloads stats::cmdscale turns on the liosting and adds slots and class attributes for which there are methods. 
#'  
#' @param d a distance structure such as that returned by 'dist' or a full symmetric matrix containing the dissimilarities
#' @param k the maximum dimension of the space which the data are to be represented in
#' @param eig indicates whether eigenvalues should be returned. Defaults to TRUE. 
#' @param ... additional parameters passed to cmdscale. See \code{\link{cmdscale}} 
#'
#' @return Object of class 'cmdscalex' and 'cmdscale' extending \code{\link{cmdscale}}. This wrapper always returns the results of cmdscale as a list, adds column labels to the $points and adds extra elements (conf=points, delta=d, confdist=dist(conf), dhat=d) and the call to the list, and assigns S3 class 'cmdscalex' and 'cmdscale'.
#'
#' @importFrom stats cmdscale as.dist dist
#'
#' 
#' @export
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#' res<-cmdscale(dis)
cmdscale <- function(d, k=2, eig=FALSE,...)
{
    if(is.data.frame(d)) d <- as.matrix(d)
    out <- stats::cmdscale(d, k=k, eig=eig, list.=TRUE,...)
    colnames(out$points) <- paste("D",1:k,sep="")
    out$conf <- out$points
    out$confdist <- stats::dist(out$conf)
    out$delta <- stats::as.dist(d)
    out$dhat <- stats::as.dist(d)
    out$call <- match.call()
    class(out) <- c("cmdscalex","cmdscale")
    out
 }

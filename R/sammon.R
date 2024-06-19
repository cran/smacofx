#'Wrapper to \code{sammon} for S3 class
#'
#' @details Overloads MASS::sammon and adds new slots and class attributes for which there are methods. 
#' 
#' @param d a distance structure such as that returned by 'dist' or a full symmetric matrix.  Data are assumed to be dissimilarities or relative distances, but must be positive except for self-distance.  This can contain missing values.
#' @param y An initial configuration. If NULL, \code{\link[smacofx]{cmdscale}} is used to provide the classical solution. (If there are missing values in 'd', an initial configuration must be provided.)  This must not have duplicates.
#' @param k The dimension of the configuration
#' @param ... Additional parameters passed to \code{sammon}, see \code{\link{sammon}}  
#'
#' @return Object of class 'sammonx' inheriting from \code{\link{sammon}}. This wrapper adds an extra slot to the list with the call, adds column labels to the $points, adds slots conf=points, delta=d, dhat=normalized dissimilarities, confdist=distance between points in conf, stress.m=stress, stress=sqrt(stress.m) and assigns S3 classes 'sammonx', 'sammon' and 'cmdscalex'.
#'
#' @importFrom MASS sammon
#' @importFrom stats as.dist dist 
#' 
#' @export
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#' res<-sammon(dis)
sammon <- function(d,y=NULL,k=2,...)
{
     if(is.data.frame(d)) d <- as.matrix(d)
     if(is.null(y)) y <- smacofx::cmdscale(d,k,eig=TRUE)$conf
     out <- MASS::sammon(d,y=as.matrix(y),k=k,...)
     colnames(out$points) <- paste("D",1:k,sep="")
     out$conf <- out$points
     out$delta <- stats::as.dist(d)
     out$dhat <- stats::as.dist(d/enorm(d))
     out$confdist <- stats::dist(out$points)
     out$stress.m <- out$stress
     out$stress <- sqrt(out$stress.m)
     out$call <- match.call()
     class(out) <- c("sammonx","sammon","cmdscalex")
     out
 }

#' Curvilinear Component Analysis (CLCA)
#'
#'
#' A wrapper to run curvilinear component analysis via \code{\link[ProjectionBasedClustering]{CCA}} and returning a 'smacofP' object. Note this functionality is rather rudimentary.   
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances. 
#' @param Epochs Scalar; gives the number of passes through the data.
#' @param lambda0 the boundary/neighbourhood parameter(s) (called lambda_y in the original paper). It is supposed to be a numeric scalar. It defaults to the 90\% quantile of delta. 
#' @param alpha0 (scalar) initial step size, 0.5 by default
#' @param weightmat not used
#' @param init starting configuration, not used
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration; not used
#' @param itmax maximum number of iterations. Not used.
#' @param verbose should iteration output be printed; not used
#' @param method Distance calculation; currently not used.
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration
#'
#' @return a 'smacofP' object. It is a list with the components
#' \itemize{
#' \item delta: Observed, untransformed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Configuration dissimilarities
#' \item conf: Matrix of fitted configuration
#' \item stress: Default stress (stress-1; sqrt of explicitly normalized stress)
#' \item spp: Stress per point 
#' \item ndim: Number of dimensions
#' \item model: Name of model
#' \item niter: Number of iterations (training length)
#' \item nobj: Number of objects
#' \item type: Type of MDS model. Only ratio here.
#' \item weightmat: weighting matrix as supplied
#' \item stress.m: Default stress (stress-1^2)
#' \item tweightmat: transformed weighting matrix; it is weightmat here.
#'}
#'
#'
#' @details
#' This implements CCA as in Demartines & Herault (1997). A different take on the ideas of curvilinear compomnent analysis is available in the experimental functions \code{\link{spmds}} and \code{\link{spmds}}.  
#' 
#'  
#' @importFrom stats dist as.dist quantile
#' @importFrom ProjectionBasedClustering CCA
#' @importFrom vegan isomapdist 
#' @export
#' 
#' @examples
#' dis<-smacof::morse
#' res<-clca(dis,lambda0=0.4)
#' res
#' summary(res)
#' plot(res)
#' 
clca <- function (delta, Epochs=20,  alpha0 = 0.5, lambda0, ndim = 2, weightmat=1-diag(nrow(delta)), init= NULL, acc=1e-06, itmax=10000, verbose=0, method = "euclidean", principal=FALSE)
{
    cc <- match.call()
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    type <- "ratio"
    labos <- rownames(delta)
    deltaorig <- delta
    if(missing(lambda0)) lambda0 <- 3*max(stats::sd(as.matrix(delta)))
    potency_curve <-  function(v0, vn, l) return(v0 * (vn/v0)^((0:(l - 1))/(l - 1)))
    n <- nrow (delta)
    p <- ndim
    lambo <- potency_curve(lambda0, 0.01, Epochs*n)
    if(verbose>0) cat("Minimizing",type,"CLCA Stress in", Epochs,"epochs with lambda0=",lambda0,"and alpha0=", alpha0,"\n")
    deltaold <- delta
    delta <- delta/enorm(delta,weightmat)
    tmp <- ProjectionBasedClustering::CCA(DataOrDistances=delta,Epochs=Epochs,alpha0=alpha0,lambda0=lambda0,OutputDimension=ndim,PlotIt=FALSE,method=method)
    out <- list()
    xnew <- tmp$ProjectedPoints #are already normalized somehow
    #xnew <- xnew/enorm(xnew) 
    ## relabeling as they were removed in the optimal scaling
    attr(xnew,"dimnames")[[1]] <- labos
    attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
    deltaorig <- stats::as.dist(deltaorig)
    deltaold <- stats::as.dist(deltaold)
    delta <- stats::as.dist(delta)
    dout <- stats::dist(xnew)
    weightmato <- stats::as.dist(weightmat)
    weightmat <- stats::as.dist(weightmat)
    #weightmat[dout>utils::tail(lambo,1)] <- 0 #Since we go through more than one lambda here, juts let it be 
    spoint <- spp(delta, dout, weightmat) 
    resmat<-spoint$resmat
    rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])
    spp <- spoint$spp
    #spp <- colMeans(resmat)
     if (principal) {
        xnew_svd <- svd(xnew)
        xnew <- xnew %*% xnew_svd$v
     }
    snew <- tmp$Error
    itel  <- Epochs*n
    if(verbose>1) cat("*** Stress:",snew, "; Stress 1 (default reported):",sqrt(snew), "\n")  
    out <- list(delta=deltaorig, dhat=delta, confdist=dout, iord=sort(delta), conf = xnew, stress=sqrt(snew), spp=spp,  ndim=p, weightmat=weightmato, resmat=resmat, rss=rss, init=init, model="CLCA", niter = itel, nobj = dim(xnew)[1], type = type, call=cc, stress.m=snew, alpha = NA, sigma = snew, tdelta=deltaold, lambo=lambo, tweightmat=weightmat)
    out$parameters  <- c(lambda0=lambda0,alpha0=alpha0)
    out$theta <- out$pars <- out$parameters
    class(out) <- c("smacofP","smacofB","smacof")
    out
   }

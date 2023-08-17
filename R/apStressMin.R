#' Approximate Power Stress SMACOF
#'
#' An implementation to minimize power stress by majorization with ratio or interval optimal scaling. Usually more accurate but slower than 'powerStressFast'. Uses a repeat loop.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1 
#' @param type what type of MDS to fit. Only "ratio" currently.
#' @param weightmat a binary matrix of finite nonegative weights.
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration
#'
#' @return a 'smacofP' object (inheriting from 'smacofB', see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed, untransformed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Configuration dissimilarities
#' \item conf: Matrix of fitted configuration
#' \item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
#' \item spp: Stress per point 
#' \item ndim: Number of dimensions
#' \item model: Name of smacof model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of MDS model 
#' \item weightmat: weighting matrix as supplied 
#' \item stress.m: Default stress (stress-1^2)
#' \item tweightmat: transformed weightingmatrix (here weightmat^nu)
#' }
#'
#' @section Note:
#' Internally we calculate the approximation parameters upsilon=nu+2*lambda*(1-(1/kappa)) and tau=lambda/kappa. They are not output. 
#'
#' @importFrom stats dist as.dist
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-apStressMin(as.matrix(dis),kappa=2,lambda=1.5,itmax=1000)
#' res
#' summary(res)
#' plot(res)
#' plot(res,"Shepard")
#' plot(res,"transplot")
#' 
#' @export
apStressMin <- function (delta, kappa=1, lambda=1, nu=1, type="ratio", weightmat= 1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    #TODO add optional arguments tau=lambda/kappa, upsilon=nu+2*lambda*(1-(1/kappa))
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(length(setdiff(unique(unlist(as.vector(weightmat))),c(0,1)))>0) stop("For approximated power stress, only binary weight matrices are allowed.")
    if(is.null(init)) init <- "torgerson"
    if(inherits(weightmat,"dist") || is.data.frame(weightmat)) weightmat <- as.matrix(weightmat)
    if(!isSymmetric(weightmat)) stop("weightmat is not symmetric.\n")
    type <- match.arg(type, c("ratio"),several.ok = FALSE) 
    trans <- "ratio"
    typo <- "none"
    if(verbose>0) cat("Minimizing approximate",type,"power-stress with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    n <- nrow (delta)
    upsilon <- nu+2*lambda*(1-(1/kappa))
    tau <- lambda/kappa
    p <- ndim
    if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
    deltaorig <- delta
    tdelta <- delta^tau
    combwght <- weightmat*(delta^upsilon)
    out <- smacof::smacofSym(tdelta,type=type,weightmat=combwght,itmax=itmax,verbose=verbose,principal=principal,init=init,ndim=ndim,eps=acc)
    tweightmat <- out$weightmat
    #stressen <- sum(weightmat*(doute-delta)^2)
    out$delta <- stats::as.dist(deltaorig)
    out$tweightmat <-stats::as.dist(tweightmat)
    out$weightmat <- stats::as.dist(weightmat)
    out$model <- "Approx. Power-Stress SMACOF"
    out$call <- match.call()
    out$stress.m <- out$stress^2
    out$tdelta <- stats::as.dist(tdelta)
    #TODO: In approx power stress if we do the kappa or it will give us. Use parameters only for print and pars for the rest  
    out$parameters <- c(kappa=kappa,lambda=lambda,nu=nu)
    out$pars <- c(kappa=kappa,lambda=lambda,nu=nu)#,upsilon=upsilon,tau=tau)
    out$theta <- c(kappa=kappa,lambda=lambda,nu=nu)
    #out$tweightmat <- weightmat
    class(out) <- c("smacofP","smacofB","smacof")
    out
}

#' @rdname apStressMin
#' @export
apowerstressMin <- apStressMin

#' @rdname apStressMin
#' @export
apostmds <- apStressMin

#' @rdname apStressMin
#' @export
apstressMin <- apStressMin

#' @rdname apStressMin
#' @export
apstressmds <- apStressMin


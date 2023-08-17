
#' Restricted Power Stress SMACOF
#'
#' An implementation to minimize restricted power stress by majorization with ratio or interval optimal scaling. Restricted means that the same power is used for both dissimilarities and fitted distances. Uses a repeat loop.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param expo power of the transformation of the fitted distances and dissimilarities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1
#' @param type what type of MDS to fit. One of "ratio" or "interval". Default is "ratio".
#' @param weightmat a matrix of finite weights
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
#' \item tweightmat: transformed weighthing matrix (here weightmat^nu)
#' \item parameters, pars, theta: The parameter vector of the explicit transformations
#' }
#'
#' 
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-rpowerStressMin(as.matrix(dis),expo=1.7,itmax=1000)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
rpowerStressMin <- function (delta, expo=1, nu=1,  type="ratio", weightmat, init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    if(missing(weightmat)) weightmat <- 1-diag(nrow(delta))
    out <- powerStressMin(delta=delta, kappa=expo, lambda=expo, nu=nu,  type=type, weightmat=weightmat, init=init, ndim = ndim, acc= acc, itmax = itmax, verbose = verbose, principal=principal)
    return(out)
}

#' @rdname rpowerStressMin
#' @export
rpowerstressMin <- rpowerStressMin

#' @rdname rpowerStressMin
#' @export
rpostmds <- rpowerStressMin

#' @rdname rpowerStressMin
#' @export
rpstressMin <- rpowerStressMin

#' @rdname rpowerStressMin
#' @export
rpStressMin <- rpowerStressMin

#' @rdname rpowerStressMin
#' @export
rpstressmds <- rpowerStressMin

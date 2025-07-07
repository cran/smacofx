
#' Power stress minimization by NEWUOA (nloptr)
#'
#' An implementation to minimize power stress by a derivative-free trust region optimization algorithm (NEWUOA). Much faster than majorizing as used in powerStressMin but perhaps less accurate. 
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1 
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc  The smallest value of the trust region radius that is allowed. If not defined, then 1e-6 will be used.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#'
#' @return a 'smacofP' object (inheriting from 'smacofB', see \code{\link[smacof]{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities, not normalized
#' \item obsdiss: Observed dissimilarities, normalized 
#' \item confdist: Configuration dissimilarities, NOT normalized 
#' \item conf: Matrix of fitted configuration, NOT normalized
#' \item stress: Default stress (stress 1, square root of the explicitly normalized stress on the normalized, transformed dissimilarities)  
#' \item spp: Stress per point (based on stress.en) 
#' \item ndim: Number of dimensions
#' \item model: Name of smacof model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of MDS model
#' }
#' and some additional components
#' \itemize{
#' \item gamma: Empty
#' \item stress.m: default stress for the COPS and STOP. Defaults to the explicitly normalized stress on the normalized, transformed dissimilarities
#' \item stress.en: explicitly stress on the normalized, transformed dissimilarities and normalized transformed distances
#' \item deltaorig: observed, untransformed dissimilarities
#' \item weightmat: weighting matrix 
#'}
#'
#' @importFrom stats dist as.dist
#' @importFrom minqa newuoa
#' 
#' @seealso \code{\link[smacof]{smacofSym}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-powerStressFast(as.matrix(dis),kappa=2,lambda=1.5)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
powerStressFast <- function (delta, kappa=1, lambda=1, nu=1, weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc = 1e-6, itmax = 10000, verbose = FALSE)
{
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(verbose>0) cat("Minimizing powerstress by NEWUOA with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")
    type <- "ratio"
    r <- kappa/2
    p <- ndim
    deltaorig <- delta
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 1 #new
    deltaold <- delta
    delta <- delta / enorm (delta, weightmat) #sum=1
    xold <- init
    if(is.null(init)) xold <- smacof::torgerson (delta, p = ndim)
    xold <- xold/enorm(xold) 
    stressf <- function(x,delta,r,ndim,weightmat)
           {
             if(!is.matrix(x)) x <- matrix(x,ncol=ndim)
             delta <- delta/enorm(delta,weightmat)
             x <- x/enorm(x)
             #adapted from powerStressMin 
             dnew <- sqdist (x)
             rnew <- sum (weightmat * delta * mkPower (dnew, r))
             nnew <- sum (weightmat * mkPower (dnew,  2*r))
             anew <- rnew / nnew
             snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
             snew
           }
     suppressWarnings(optimized <- minqa::newuoa(xold,function(par) stressf(par,delta=delta,r=r,ndim=ndim,weightmat=weightmat),control=list(maxfun=itmax,rhoend=acc,iprint=verbose-2)))
      xnew <- matrix(optimized$par,ncol=ndim)
      xnew <- xnew/enorm(xnew)
             #adapted from powerStressMin 
      dnew <- sqdist (xnew)
      rnew <- sum (weightmat * delta * mkPower (dnew, r))
      nnew <- sum (weightmat * mkPower (dnew,  2*r))
      anew <- rnew / nnew
      stress <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
     #xnew <- optimized$par
      attr(xnew,"dimnames")[[1]] <- rownames(delta)
      itel <- optimized$feval
      attr(xnew,"dimnames")[[2]] <- paste("D",1:ndim,sep="")
      doutm <- (2*sqrt(sqdist(xnew)))^kappa  #fitted powered euclidean distance but times two
      deltam <- delta
      deltaorigm <- deltaorig
      deltaoldm <- deltaold
      delta <- stats::as.dist(delta)
      deltaorig <- stats::as.dist(deltaorig)
      deltaold <- stats::as.dist(deltaold)
      doute <- doutm/enorm(doutm)
      doute <- stats::as.dist(doute)
      dout <- stats::as.dist(doutm)
      resmat <- as.matrix(delta - doute)^2
      spp <- colMeans(resmat)
      weightmatm <-weightmat
      weightmat <- stats::as.dist(weightmatm)
      stressen <- sum(weightmat*(doute-delta)^2) #raw stress on the normalized proximities and normalized distances 
      if(verbose>1) cat("*** stress (both normalized - for COPS/STOPS):",stress,"; stress 1 (both normalized - default reported):",sqrt(stress),"; stress manual (for debug only):",stressen,"; from optimization: ",optimized$fval,"\n")   
    out <- list(delta=deltaold, obsdiss=delta, confdist=dout, conf = xnew, pars=c(kappa,lambda,nu), niter = itel, stress=stress, spp=spp, ndim=p, model="power-stress MDS (NEWUOA)", call=match.call(), nobj = dim(xnew)[1], type = type, gamma = NA, stress.m=sqrt(stress), stress.en=stressen, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat)
    class(out) <- c("smacofP","smacofB","smacof")
    out
}

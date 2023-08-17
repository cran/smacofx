#' Sparsified (POST-) Multidimensional Scaling (SPMDS or SMDS) either as self-organizing or not
#'
#' An implementation of a sparsified version of (POST-)MDS by pseudo-majorization with ratio, interval and ordinal optimal scaling for dissimilarities and optional power transformations. This is inspired by curvilinear component analysis but works differently: It finds an initial weightmatrix where w_ij(X^0)=0 if d_ij(X^0)>tau and fits a POST-MDS with these weights. Then in each successive iteration step, the weightmat is recalculated so that w_ij(X^(n+1))=0 if d_ij(X^(n+1))>tau. 
#'
#' There is a wrapper 'smds' where the exponents are 1, which is standard SMDS but extend to allow optimal scaling. The neighborhood parameter tau is kept fixed in 'spmds' and 'smds'. The functions 'so_spmds' and 'so_smds' implement a self-organising principle, where the SMDS is repeatedly fitted for a decreasing sequence of taus.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param lambda exponent of the power transformation of the dissimilarities; defaults to 1, which is also the setup of 'smds'
#' @param kappa exponent of the power transformation of the fitted distances; defaults to 1, which is also the setup of 'smds'.
#' @param nu exponent of the power of the weighting matrix; defaults to 1 which is also the setup for 'smds'. 
#' @param tau the boundary/neighbourhood parameter(s) (called lambda in the original paper). For 'spmds' and 'smds' it is supposed to be a numeric scalar (if a sequence is supplied the maximum is taken as tau) and all the transformed fitted distances exceeding tau are set to 0 via the weightmat (assignment can change between iterations). It defaults to the 90\% quantile of delta. For 'so_spmds' tau is supposed to be either a user supplied decreasing sequence of taus or if a scalar the maximum tau from which a decreasing sequence of taus is generated automatically as 'seq(from=tau,to=tau/epochs,length.out=epochs)' and then used in sequence.
#' @param type what type of MDS to fit. Currently one of "ratio", "interval" or "ordinal". Default is "ratio".
#' @param ties the handling of ties for ordinal (nonmetric) MDS. Possible are "primary" (default), "secondary" or "tertiary".
#' @param weightmat a matrix of finite weights. 
#' @param init starting configuration. If NULL (default) we fit a full rstress model.
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration
#' @param epochs for 'so_spmds' and tau being scalar, it gives the number of passes through the data. The sequence of taus created is 'seq(tau,tau/epochs,length.out=epochs)'. If tau is of length >1, this argument is ignored.
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
#' \item tweightmat: transformed weighting matrix; it is weightmat but containing all the 0s for the distances set to 0. 
#'}
#'
#'
#' @details
#' The solution is found by "quasi-majorization", which means that the majorization is only real majorization once the weightmat no longer changes. This typically happens after a few iterations. Due to that it can be that in the beginning the stress may not decrease monotonically and that there's a chance it might never. 
#' 
#' If tau is too small it may happen that all distances for one i to all j are zero and then there will be an error, so make sure to set a larger tau.
#'
#' In the standard functions 'spmds' and 'smds' we keep tau fixed throughout. This means that if tau is large enough, then the result is the same as the corresponding MDS. In the orginal publication the idea was that of a self-organizing map which decreased tau over epochs (i.e., passes through the data). This can be achieved with our function 'so_spmds' 'so_smds' which creates a vector of decreasing tau values, calls the function 'spmds' with the first tau, then supplies the optimal configuration obtained as the init for the next call with the next tau and so on. 
#'
#' 
#' @importFrom stats dist as.dist quantile
#' @importFrom smacof transform transPrep
#' 
#' @examples
#' dis<-smacof::morse
#' res<-spmds(dis,type="interval",kappa=2,lambda=2,tau=0.35,itmax=100) #use higher itmax
#' res2<-smds(dis,type="interval",tau=0.35,itmax=500) #use higher itmax
#' res
#' res2
#' summary(res)
#' oldpar<-par(mfrow=c(1,2))
#' plot(res)
#' plot(res2)
#' par(oldpar)
#'
#' ##which d_{ij}(X) exceeded tau at convergence (i.e., have been set to 0)?
#' res$tweightmat
#' res2$tweightmat
#'
#' \donttest{
#' ## Self-organizing map style (as in the clca publication)
#' #run the som-style (p)smds 
#' sommod1<-so_spmds(dis,tau=0.3,kappa=0.5,lambda=2,epochs=10,verbose=1)
#' sommod2<-so_smds(dis,tau=0.3,epochs=10,verbose=1)
#' sommod1
#' sommod2
#' }
#' 
#' @export
spmds <- function (delta, lambda=1, kappa=1, nu=1, tau, type="ratio", ties="primary", weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    if(inherits(weightmat,"dist") || is.data.frame(weightmat)) weightmat <- as.matrix(weightmat)
    if(!isSymmetric(weightmat)) stop("weightmat is not symmetric.\n")
    r <- kappa/2
    if(length(tau)>1)
    {
        warning("Supplied tau is of length >1. The max(tau) was used as tau.")
        tau <- max(tau)
    }
    if(tau<=0) stop("tau must be positive.")
    ## -- Setup for MDS type
    if(missing(type)) type <- "ratio"
    type <- match.arg(type, c("ratio", "interval", "ordinal"),several.ok = FALSE)
    #if(type =="ordinal") lambda <- 1 #We dont allow powers for dissimilarities in nonmetric MDS
    #    "mspline"), several.ok = FALSE)
    trans <- type
    typo <- type
    if (trans=="ratio"){
    trans <- "none"
    }
    else if (trans=="ordinal" & ties=="primary"){
    trans <- "ordinalp"
    typo <- "ordinal (primary)"
   } else if(trans=="ordinal" & ties=="secondary"){
    trans <- "ordinals"
    typo <- "ordinal (secondary)"
  } else if(trans=="ordinal" & ties=="tertiary"){
    trans <- "ordinalt"
    typo <- "ordinal (tertiary)"
  #} else if(trans=="spline"){
  #  trans <- "mspline"
  }
    if(verbose>0) cat(paste("Fitting",type,"spmds with lambda=",lambda, "kappa=",kappa,"nu=",nu, "and tau=",tau,"\n"))
    n <- nrow (delta)
    normi <- 0.5
    ##normi <- n #if normi=n we can use the iord structure in plot.smacofP
    ## but the problem is we don't get the correct stress then anymore.
    p <- ndim
    if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
    if(is.null(rownames(delta))) rownames(delta) <- 1:n 
    labos <- rownames(delta) #labels
    deltaorig <- delta
    delta <- delta^lambda
    weightmato <- weightmat
    weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 0
    delta <- delta / enorm (delta, weightmat)
    if(missing(tau)) tau <- stats::quantile(delta,0.9)
    disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)#spline.intKnots = spline.intKnots, spline.degree = spline.degree) #FIXME: only works with dist() style object 
    ## Add an intercept to the spline base transformation
                                        #if (trans == "mspline") disobj$base <- cbind(rep(1, nrow(disobj$base)), disobj$base)
    #delta <- delta / enorm (delta, weightmat)
    deltaold <- delta
    itel <- 1
    ##Starting Configs
    xold  <- init
    if(is.null(init)) xold <- smacofx::rStressMin (delta,r=kappa/2,type=type,ties=ties,weightmat=weightmat,ndim=ndim,init=init,itmax=itmax,principal=principal)$conf
    xstart <- xold
    xold <- xold / enorm (xold) 
    nn <- diag (n)
    dold <- sqdist (xold)
    weightmat[sqrt(dold)>tau] <- 0 ##CCA penalty
    ##first optimal scaling
    eold <- as.dist(sqrt(dold))
    dhat <- smacof::transform(eold, disobj, w = as.dist(weightmat), normq = normi)
    dhatt <- dhat$res #I need the structure here to reconstruct the delta; alternatively turn all into vectors? - checked how they do it in smacof
    dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
    #FIXME: labels
    delta <- as.matrix(dhatd)
    rold <- sum (weightmat * delta * mkPower (dold, r))
    nold <- sum (weightmat * mkPower (dold, 2 * r))
    aold <- rold / nold
    sold <- 1 - 2 * aold * rold + (aold ^ 2) * nold
    ## Optimizing
    repeat {
      if(tau<=min(sqrt(dold[lower.tri(dold)]))) stop("Current tau is lower than the smallest fitted distance (so all distances are set to 0). Increase tau.")
      p1 <- mkPower (dold, r - 1)
      p2 <- mkPower (dold, (2 * r) - 1)
 
      by <- mkBmat (weightmat * delta * p1)
      cy <- mkBmat (weightmat * p2)
      ga <- 2 * sum (weightmat * p2)
      be <- (2 * r - 1) * (2 ^ r) * sum (weightmat * delta)
      de <- (4 * r - 1) * (4 ^ r) * sum (weightmat)
      if (r >= 0.5) {
        my <- by - aold * (cy - de * nn)
      }
      if (r < 0.5) {
        my <- (by - be * nn) - aold * (cy - ga * nn)
      }
      xnew <- my %*% xold
      xnew <- xnew / enorm (xnew)
      dnew <- sqdist (xnew)
      ## TODO: What do we need? Two possibilities
      ### Should we always set the 0 freshly, so it is possible that a w_{ij} can change from 0 to >0 again? then we need next line     
      weightmat <- weightmato #new
      ## or should we never change the 0 back once it was found? I'm sure that then the algorithm is majorizing this objective; it also coincides with the above if the d_{ij] are monotonically decreasing but I don't think that is what they do.
      ## test this 
      weightmat[!is.finite(weightmat)] <- 0
      weightmat[sqrt(dnew)>tau] <- 0
      ##optimal scaling
      e <- as.dist(sqrt(dnew)) #I need the dist(x) here for interval
      dhat2 <- smacof::transform(e, disobj, w = as.dist(weightmat), normq = normi)  ## dhat update
      dhatt <- dhat2$res 
      dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
      delta <- as.matrix(dhatd)
      #delta <- as.matrix(dhatt) #In cops this is <<- because we need to change it outside of copsf() but here we don't need that 
      rnew <- sum (weightmat * delta * mkPower (dnew, r))
      nnew <- sum (weightmat * mkPower (dnew, 2 * r))
      anew <- rnew / nnew
      snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
      if(is.na(snew)) #if there are issues with the values
          {
              snew <- sold
              dnew <- dold
              anew <- aold
              xnew <- xold
          }   
      if (verbose>2) {
        cat (
          formatC (itel, width = 4, format = "d"),
          formatC (
            sold,
            digits = 10,
            width = 13,
            format = "f"
          ),
          formatC (
            snew,
            digits = 10,
            width = 13,
            format = "f"
          ),
          "\n"
        )
      }
      if ((itel == itmax) || (abs(sold - snew) < acc)) #new
        break ()
      itel <- itel + 1
      xold <- xnew
      dold <- dnew
      sold <- snew
      aold <- anew
    }
    xnew <- xnew/enorm(xnew)
    ## relabeling as they were removed in the optimal scaling
    rownames(delta) <- labos
    attr(xnew,"dimnames")[[1]] <- rownames(delta)
    attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
    doutm <- mkPower(sqdist(xnew),r)
    deltam <- delta
    #delta <- structure(delta, Size = n, call = quote(as.dist.default(m=b)),
    #                   class = "dist", Diag = FALSE, Upper = FALSE)
    delta <- stats::as.dist(delta)
    deltaorig <- stats::as.dist(deltaorig)
    deltaold <- stats::as.dist(deltaold)
    #doute <- doutm/enorm(doutm) #this is an issue here!
    #doute <- stats::as.dist(doute)
    dout <- stats::as.dist(doutm)
    weightmatm <-weightmat
    #resmat <- weightmatm*as.matrix((delta - doute)^2) #Old version 
    #resmat <- weightmatm*as.matrix((deltam - doutm)^2)
    weightmat <- stats::as.dist(weightmatm)
    #spp <- colMeans(resmat)
    spoint <- spp(delta, dout, weightmat)
    resmat<-spoint$resmat
    rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])
    spp <- spoint$spp
    #spp <- colMeans(resmat)
    if (principal) {
        xnew_svd <- svd(xnew)
        xnew <- xnew %*% xnew_svd$v
    }
    #stressen <- sum(weightmat*(doute-delta)^2)
    if(verbose>1) cat("*** Stress:",snew, "; Stress-1 (default reported):",sqrt(snew),"\n")
    #delta is input delta, tdelta is input delta with explicit transformation and normalized, dhat is dhats 
    out <- list(delta=deltaorig, dhat=delta, confdist=dout, iord=dhat2$iord.prim, conf = xnew, stress=sqrt(snew), spp=spp,  ndim=p, weightmat=weightmato, resmat=resmat, rss=rss, init=xstart, model="power SMDS", niter = itel, nobj = dim(xnew)[1], type = type, call=match.call(), stress.m=snew, alpha = anew, sigma = snew, tdelta=deltaold, parameters=c(kappa=kappa,lambda=lambda,nu=nu,tau=tau), pars=c(kappa=kappa,lambda=lambda,nu=nu,tau=tau), theta=c(kappa=kappa,lambda=lambda,nu=nu,tau=tau),tweightmat=weightmat)
    class(out) <- c("smacofP","smacofB","smacof")
    out
  }


#' @rdname spmds
#' @export
smds <- function(delta, tau=stats::quantile(delta,0.9), type="ratio", ties="primary", weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    cc <- match.call()
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    out <- spmds(delta=delta, lambda=1, kappa=1, nu=1, tau=tau, type=type, ties=ties, weightmat=weightmat, init=init, ndim=ndim, acc=acc, itmax=itmax, verbose=verbose, principal=principal)
    out$model <- "SMDS"
    out$call <- cc
    out$parameters <- out$theta <- out$pars  <- c(tau=tau)
    out
}

#' @rdname spmds
#' @export
so_spmds <- function(delta, kappa=1, lambda=1, nu=1, tau=max(delta), epochs=10, type="ratio", ties="primary", weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    cc <- match.call()
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    if(length(tau)<2)
       {
         taumax <- tau
         taumin <- tau/epochs
         taus <- seq(taumax,taumin,length.out=epochs)
       } else taus <- tau
    if(any(diff(taus)>0)) taus <- sort(taus,decreasing=TRUE)
    finconf <- init
    for(i in 1:length(taus))
    {
      if(verbose>0) cat(paste0("Epoch ",i,": tau=",taus[i],"\n"))  
      tmp<-spmds(delta=delta, lambda=lambda, kappa=kappa, nu=nu, tau=taus[i], type=type, ties=ties, weightmat=weightmat, init=finconf, ndim=ndim, verbose=verbose-1, acc=acc, itmax=itmax, principal=principal)
      finconf<-tmp$conf
      finmod<-tmp
    }
    finmod$call  <- cc
    finmod$model  <- "SO-SPMDS"
    return(finmod)
}

#' @rdname spmds
#' @export
so_smds <- function(delta, tau=max(delta), epochs=10, type="ratio", ties="primary", weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    cc <- match.call()
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    if(length(tau)<2)
       {
         taumax <- tau
         taumin <- tau/epochs
         taus <- seq(taumax,taumin,length.out=epochs)
       } else taus <- tau
    if(any(diff(taus)>0)) taus <- sort(taus,decreasing=TRUE)
    finconf <- init
    for(i in 1:length(taus))
    {
      if(verbose>0) cat(paste0("Epoch ",i,": tau=",taus[i],"\n"))  
      tmp<-smds(delta=delta, tau=taus[i], type=type, ties=ties, weightmat=weightmat, init=finconf, ndim=ndim, verbose=verbose,  acc=acc, itmax=itmax, principal=principal)
      finconf<-tmp$conf
      finmod<-tmp
    }
    finmod$call  <- cc
    finmod$model  <- "SO-SMDS"
    return(finmod)
    }

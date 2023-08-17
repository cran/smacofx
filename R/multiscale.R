## #' Multiscale SMACOF
## #'
## #' An implementation for maximum likleihood MDS aka multiscale that minimizes the multiscale stress by majorization with ratio, interval and ordinal optimal scaling. Uses a repeat loop.
## #' 
## #' @param delta dist object or a symmetric, numeric data.frame or matrix of distances. Warning: these will get transformed to the log scale, so make sure that log(delta)>=0.    
## #' @param type what optimal scaling type of MDS to fit. Currently one of "ratio" or "interval". Default is "ratio". 
## #' @param weightmat a matrix of finite weights
## #' @param init starting configuration
## #' @param ndim dimension of the configuration; defaults to 2
## #' @param acc numeric accuracy of the iteration. Default is 1e-6.
## #' @param itmax maximum number of iterations. Default is 10000.
## #' @param verbose should iteration output be printed; if > 1 then yes
## #' @param kappa As this is not exactly multiscale but an r-stress approximation, we have multiscale only for kappa->0. This argument can therefore be used to make the approximation more accurate by making it smaller. Default is 0.1.
## #' 
## #' @return a smacofP object (inheriting form smacofB, see \code{\link{smacofSym}}). It is a list with the components
## #' \itemize{
## #' \item delta: Observed dissimilarities on log scale 
## #' \item obsdiss: Observed dissimilarities (dhats), optimally scaled and normalized 
## #' \item confdist: Configuration dissimilarities, NOT normalized 
## #' \item conf: Matrix of fitted configuration, NOT normalized
## #' \item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
## #' \item spp: Stress per point (based on stress.en) 
## #' \item ndim: Number of dimensions
## #' \item model: Name of smacof model
## #' \item niter: Number of iterations
## #' \item nobj: Number of objects
## #' \item type: Type of MDS model
## #' }
## #' and some additional components
## #' \itemize{
## #' \item stress.m: default stress for the COPS and STOP defaults to the explicitly normalized stress on the normalized, transformed dissimilarities
## #' \item deltaorig: observed, untransformed dissimilarities (as inputted)
## #' \item weightmat: weighting matrix 
## #'}
## #'
## #' @section
## #' Warning: The input delta will internally get transformed to the log scale, so make sure that log(delta)>=0 otherwise it throws an error. It is often a good idea to use 1+delta in this case. 
## #' 
## #'
## #' @importFrom stats dist as.dist
## #' 
## #' @seealso \code{\link{rStressMin}}
## #' 
## #' @examples
## #' dis<-smacof::kinshipdelta
## #' res<-multiscale(as.matrix(dis),type="interval",itmax=1000)
## #' res
## #' summary(res)
## #' plot(res)
## #' 
## #' @export
## multiscaleOLD <- function (delta, type=c("ratio","interval"), weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE,kappa=0.1) {
##     if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
##     if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
##     r <- kappa/2
##     ## -- Setup for MDS type
##     if(missing(type)) type <- "ratio"
##     trans <- type
##     typo <- type
##     if (trans=="ratio"){
##     trans <- "none"
##     }
##    # else if (trans=="ordinal" & ties=="primary"){
##    # trans <- "ordinalp"
##    # typo <- "ordinal (primary)"
##    #} else if(trans=="ordinal" & ties=="secondary"){
##    # trans <- "ordinals"
##    # typo <- "ordinal (secondary)"
##   #} else if(trans=="ordinal" & ties=="tertiary"){
##   #  trans <- "ordinalt"
##   #  typo <- "ordinal (tertiary)"
##   #} else if(trans=="spline"){
##   #  trans <- "mspline"
##   #}
##     if(verbose>0) cat(paste("Minimizing",type,"rStress with r=",r,"\n"))    
##     p <- ndim
##     labos <- rownames(delta) #labels
    
##     deltaorig <- delta
##     delta <- log(delta) #The multiscale dissimilarities 
##     diag(delta) <- 0 #making main diagonal 0.
##     if(any(delta)<0) stop("Nonnegative dhats after log transformation. Input different dissimilarities.\n")
##     #weightmato <- weightmat
##     #weightmat <- weightmat^nu
##     weightmat[!is.finite(weightmat)] <- 1
##     deltaold <- delta
##     delta <- delta / enorm (delta, weightmat)
##     disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)#spline.intKnots = spline.intKnots, spline.degree = spline.degree) #FIXME: only works with dist() style object 
##     ## Add an intercept to the spline base transformation
##     #if (trans == "mspline") disobj$base <- cbind(rep(1, nrow(disobj$base)), disobj$base)
##     itel <- 1
##     ##Starting Configs
##     xold <- init
##     if(is.null(init)) xold <- smacof::torgerson (delta, p = p)
##     xold <- xold / enorm (xold) 
##     n <- nrow (xold)
##     nn <- diag (n)
##     dold <- sqdist (xold)
##    ##first optimal scaling
##     eold <- as.dist(sqrt(dold))
##     dhat <- smacof::transform(eold, disobj, w = as.dist(weightmat), normq = 0.5)
##     dhatt <- dhat$res #I need the structure here to reconstruct the delta; alternatively turn all into vectors? - checked how they do it in smacof
##     dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
##     #FIXME: labels
##     delta <- as.matrix(dhatd)
##     rold <- sum (weightmat * delta * mkPower (dold, r))
##     nold <- sum (weightmat * mkPower (dold, 2 * r))
##     aold <- rold / nold
##     sold <- 1 - 2 * aold * rold + (aold ^ 2) * nold
##     ## Optimizing
##     repeat {
##       p1 <- mkPower (dold, r - 1)
##       p2 <- mkPower (dold, (2 * r) - 1)
##       by <- mkBmat (weightmat * delta * p1)
##       cy <- mkBmat (weightmat * p2)
##       ga <- 2 * sum (weightmat * p2)
##       be <- (2 * r - 1) * (2 ^ r) * sum (weightmat * delta)
##       de <- (4 * r - 1) * (4 ^ r) * sum (weightmat)
##       if (r >= 0.5) {
##         my <- by - aold * (cy - de * nn)
##       }
##       if (r < 0.5) {
##         my <- (by - be * nn) - aold * (cy - ga * nn)
##       }
##       xnew <- my %*% xold
##       xnew <- xnew / enorm (xnew)
##       dnew <- sqdist (xnew)
##       ##optimal scaling
##       e <- as.dist(sqrt(dnew)) #I need the dist(x) here for interval
##       dhat <- smacof::transform(e, disobj, w = as.dist(weightmat), normq = 0.5)  ## dhat update
##       dhatt <- dhat$res #FIXME: I need the structure here to reconstruct the delta; alternatively turn all into vectors - check how they do it in smacof
##       dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
##       delta <- as.matrix(dhatd) #In cops this is <<- because we need to change it outside of copsf() but here we don't need that 
##       rnew <- sum (weightmat * delta * mkPower (dnew, r))
##       nnew <- sum (weightmat * mkPower (dnew, 2 * r))
##       anew <- rnew / nnew
##       snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
##       if(is.na(snew)) #if there are issues with the values
##           {
##               snew <- sold
##               dnew <- dold
##               anew <- aold
##               xnew <- xold
##           }   
##       if (verbose>2) {
##         cat (
##           formatC (itel, width = 4, format = "d"),
##           formatC (
##             sold,
##             digits = 10,
##             width = 13,
##             format = "f"
##           ),
##           formatC (
##             snew,
##             digits = 10,
##             width = 13,
##             format = "f"
##           ),
##           "\n"
##         )
##       }
##       if ((itel == itmax) || ((sold - snew) < acc))
##         break ()
##       itel <- itel + 1
##       xold <- xnew
##       dold <- dnew
##       sold <- snew
##       aold <- anew
##     }
##     xnew <- xnew/enorm(xnew) #TODO: Check in smacof whether I need this
##     ## relabeling as they were removed in the optimal scaling
##     rownames(delta) <- labos
##     attr(xnew,"dimnames")[[1]] <- rownames(delta)
##     attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
##     doutm <- mkPower(sqdist(xnew),r)
##     deltam <- delta
##     delta <- stats::as.dist(delta)
##     deltaorig <- stats::as.dist(deltaorig)
##     deltaold <- stats::as.dist(deltaold)
##     #doute <- doutm/enorm(doutm) #this is an issue here!
##     #doute <- stats::as.dist(doute)
##     dout <- stats::as.dist(doutm)
##     weightmatm <-weightmat
##     resmat <- weightmatm*as.matrix((deltam - doutm)^2) 
##     spp <- colMeans(resmat) 
##     weightmat <- stats::as.dist(weightmatm)
##     #stressen <- sum(weightmat*(doute-delta)^2)
##     if(verbose>1) cat("***Stress:",snew, "; Stress 1 (both normalized - default reported):",sqrt(snew),"\n")  
##     out <- list(delta=deltaold, obsdiss=delta, confdist=dout, conf = xnew, parameters=c(kappa=kappa,lambda=1), pars=c(kappa=kappa,lambda=1), theta=c(kappa=kappa,lambda=1), niter = itel, spp=spp, ndim=p, model="Multiscale SMACOF", call=match.call(), nobj = dim(xnew)[1], type = type, stress=sqrt(snew), stress.m=snew, deltaorig=as.dist(deltaorig),resmat=resmat,weightmat=weightmat, alpha = anew, sigma = snew)
##     class(out) <- c("smacofP","smacofB","smacof")
##     out
##   }

#' Multiscale SMACOF
#'
#' An implementation for maximum likelihood MDS aka multiscale that minimizes the multiscale stress by majorization with ratio and interval optimal scaling. Uses a repeat loop.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances. Warning: these will get transformed to the log scale, so make sure that log(delta)>=0.    
#' @param type what optimal scaling type of MDS to fit. Currently one of "ratio" or "interval". Default is "ratio". 
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param kappa As this is not exactly multiscale but an r-stress approximation, we have multiscale only for kappa->0. This argument can therefore be used to make the approximation more accurate by making it smaller. Default is 0.1.
#' @param principal If ‘TRUE’, principal axis transformation is applied to the final configuration
#' 
#' @return a 'smacofP' object (inheriting from 'smacofB', see \code{\link{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities
#' \item tdelta: Observed explicitly transformed (log) dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Configuration dissimilarities, NOT normalized 
#' \item conf: Matrix of fitted configuration, NOT normalized
#' \item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
#' \item spp: Stress per point
#' \item ndim: Number of dimensions
#' \item model: Name of smacof model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of MDS model
#' \item weightmat: weighting matrix 
#' \item stress.m: Default stress (stress-1^2)
#' }
#'
#' @section
#' Warning: The input delta will internally get transformed to the log scale, so make sure that log(delta)>=0 otherwise it throws an error. It is often a good idea to use 1+delta in this case. 
#' 
#'
#' @importFrom stats dist as.dist
#' 
#' @seealso \code{\link{rStressMin}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-multiscale(as.matrix(dis),type="interval",itmax=1000)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
multiscale <- function (delta, type=c("ratio","interval"), weightmat, init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, kappa=0.1, principal=FALSE) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(missing(weightmat))  weightmat <- 1-diag(nrow(delta))
    if(inherits(weightmat,"dist") || is.data.frame(weightmat)) weightmat <- as.matrix(weightmat)
    if(!isSymmetric(weightmat)) stop("weightmat is not symmetric.\n")
    r <- kappa/2
    ## -- Setup for MDS type
    if(missing(type)) type <- "ratio"
    type <- match.arg(type, c("ratio", "interval"),several.ok = FALSE)
    trans <- type
    typo <- type
    if (trans=="ratio"){
    trans <- "none"
    }
   # else if (trans=="ordinal" & ties=="primary"){
   # trans <- "ordinalp"
   # typo <- "ordinal (primary)"
   #} else if(trans=="ordinal" & ties=="secondary"){
   # trans <- "ordinals"
   # typo <- "ordinal (secondary)"
  #} else if(trans=="ordinal" & ties=="tertiary"){
  #  trans <- "ordinalt"
  #  typo <- "ordinal (tertiary)"
  #} else if(trans=="spline"){
  #  trans <- "mspline"
  #}
    if(verbose>0) cat(paste("Minimizing",type,"multiscale stress with r=",r,"\n"))    
    n <- nrow (delta)
    p <- ndim
    if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
    if(is.null(rownames(delta))) rownames(delta) <- 1:n 
    labos <- rownames(delta) #labels
    
    deltaorig <- delta
    delta <- log(delta) #The multiscale dissimilarities 
    diag(delta) <- 0 #making main diagonal 0.
    if(any(delta)<0) stop("Nonnegative delta after log transformation. Input different dissimilarities.\n")
    #weightmato <- weightmat
    #weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 0
    delta <- delta / enorm (delta, weightmat)
    disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)#spline.intKnots = spline.intKnots, spline.degree = spline.degree) #FIXME: only works with dist() style object 
    ## Add an intercept to the spline base transformation
                                        #if (trans == "mspline") disobj$base <- cbind(rep(1, nrow(disobj$base)), disobj$base)
    deltaold <- delta
    itel <- 1
    ##Starting Configs
    xold <- init
    if(is.null(init)) xold <- smacof::torgerson (delta, p = p)
    xstart <- xold
    xold <- xold / enorm (xold) 
    n <- nrow (xold)
    nn <- diag (n)
    dold <- sqdist (xold)
   ##first optimal scaling
    eold <- as.dist(sqrt(dold))
    dhat <- smacof::transform(eold, disobj, w = as.dist(weightmat), normq = 0.5)
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
      ##optimal scaling
      e <- as.dist(sqrt(dnew)) #I need the dist(x) here for interval
      dhat2 <- smacof::transform(e, disobj, w = as.dist(weightmat), normq = 0.5)  ## dhat update
      dhatt <- dhat2$res #FIXME: I need the structure here to reconstruct the delta; alternatively turn all into vectors - check how they do it in smacof
      dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
      delta <- as.matrix(dhatd) #In cops this is <<- because we need to change it outside of copsf() but here we don't need that 
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
      if ((itel == itmax) || ((sold - snew) < acc))
        break ()
      itel <- itel + 1
      xold <- xnew
      dold <- dnew
      sold <- snew
      aold <- anew
    }
    xnew <- xnew/enorm(xnew) #TODO: Check in smacof whether I need this
    ## relabeling as they were removed in the optimal scaling
    rownames(delta) <- labos
    attr(xnew,"dimnames")[[1]] <- rownames(delta)
    attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
    doutm <- mkPower(sqdist(xnew),r)
    deltam <- delta
    delta <- stats::as.dist(delta)
    deltaorig <- stats::as.dist(deltaorig)
    deltaold <- stats::as.dist(deltaold)
    #doute <- doutm/enorm(doutm) #this is an issue here!
    #doute <- stats::as.dist(doute)
    dout <- stats::as.dist(doutm)
    weightmatm <-weightmat
    #resmat <- weightmatm*as.matrix((deltam - doutm)^2) 
    #spp <- colMeans(resmat)
    weightmat <- stats::as.dist(weightmatm)
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
    if(verbose>1) cat("*** Stress:",snew, "; Stress 1 (default reported):",sqrt(snew),"\n")  
    out <- list(delta=deltaorig, dhat=delta, confdist=dout, iord=dhat2$iord.prim, conf = xnew, stress=sqrt(snew), spp=spp,  ndim=p, weightmat=weightmat, resmat=resmat, rss=rss, init=xstart, model="Multiscale SMACOF", niter = itel, nobj = dim(xnew)[1], type = type, call=match.call(), stress.m=snew, alpha = anew, sigma = snew, tdelta=deltaold, parameters=c(kappa=kappa,TDelta="log"), pars=c(kappa=kappa,TDelta="log"), theta=c(kappa=kappa,TDelta="log"),tweightmat=NULL)
    class(out) <- c("smacofP","smacofB","smacof")
    out
  }





#' Sammon Mapping SMACOF
#'
#' An implementation to minimize Sammon stress by majorization with ratio and interval optimal scaling. Uses a repeat loop.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param type what type of MDS to fit. Currently one of "ratio" and "interval". Default is "ratio".
#' @param weightmat a matrix of finite weights
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should fitting information be printed; if > 0 then yes
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration
#'
#' @return a 'smacofP' object (inheriting from smacofB, see \code{\link[smacof]{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Observed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Transformed configuration distances
#' \item conf: Matrix of fitted configuration
#' \item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
#' \item spp: Stress per point (based on stress.en) 
#' \item ndim: Number of dimensions
#' \item model: Name of smacof model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of MDS model
#' \item weightmat: weighting matrix as supplied
#' \item stress.m: default stress (stress-1^2)
#' \item tweightmat: weighting matrix atfer transformation (here weightmat/delta)
#' }
#'
#' @importFrom stats dist as.dist
#' 
#' @seealso \code{\link{rStressMin}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-sammonmap(as.matrix(dis),itmax=1000)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
sammonmap <- function (delta, type=c("ratio","interval"), weightmat, init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE,principal=FALSE) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    if(missing(weightmat))  weightmat <- 1-diag(nrow(delta))
    if(inherits(weightmat,"dist") || is.data.frame(weightmat)) weightmat <- as.matrix(weightmat)
    if(!isSymmetric(weightmat)) stop("weightmat is not symmetric.\n")
    r <- 0.5
    ## -- Setup for MDS type
    if(missing(type)) type <- "ratio"
    type <- match.arg(type, c("ratio", "interval"),several.ok = FALSE) 
    trans <- type
    typo <- type
    if (trans=="ratio"){
    trans <- "none"
    }
   ## ordinal makes no sense in sammon mapping imo
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
    if(verbose>0) cat(paste("Minimizing",type,"Sammon stress","\n"))
     n <- nrow (delta)
    p <- ndim
    if (p > (n - 1)) 
    stop("Maximum number of dimensions is n-1!")
    if(is.null(rownames(delta))) rownames(delta) <- 1:n 
    labos <- rownames(delta) #labels
   
    weightmatorig <- weightmat #back up weightmat
    deltaorig <- delta #backup delta
    weightmato <- weightmat
    ## We do the sammon weighting before anything, then normalize the delta accordingly 
    ## I checked whether we should do a first (original delta) and second (enormed delta) weighting step (we can for ratio and interval). We now use enormed delta.
    delta <- delta / enorm (delta, weightmat) #Normalize with Sammon weighted
    weightmat <- weightmato/delta #Sammon weighting #1
    weightmat[!is.finite(weightmat)] <- 0
    disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)#spline.intKnots = spline.intKnots, spline.degree = spline.degree) #FIXME: only works with dist() style object
    deltaold <- delta
    ## Add an intercept to the spline base transformation
    #if (trans == "mspline") disobj$base <- cbind(rep(1, nrow(disobj$base)), disobj$base)
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
    eold <- as.dist(mkPower(dold,r))
    ##Second samon weighting step with the Samon weighting normalized delta as w/delta * (delta-d(X)) 
    #weightmat <- weightmato/delta # Sammon weighting #2
    #weightmat[!is.finite(weightmat)] <- 0
    #we do not weight with optimally scaled dhats 
    dhat <- smacof::transform(eold, disobj, w = as.dist(weightmat), normq = 0.5) #calculate dhats. There is the sammon mapping already included via the w argument.
    dhatt <- dhat$res 
    dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
    delta <- as.matrix(dhatd)
    ##Sammon weighting #3: We don't do that because majorization isn't working anymore if we weight with optimally scaled delta 
    #weightmat <- weightmato/delta #sammon weighting with dhat and we always use weightmato/delta
    #weightmat[!is.finite(weightmat)] <- 0
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
      e <- as.dist(mkPower(dnew,r)) #I need the dist(x) here for interval
      dhat2 <- smacof::transform(e, disobj, w = as.dist(weightmat), normq = 0.5)  ## dhat update
      dhatt <- dhat2$res 
      dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
      delta <- as.matrix(dhatd) 
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
      ##Sammon mapping #4 with dynamic updating 
      ##we can't do this with interval (only two iterations then) and for ratio it makes no difference whether we have that or not so we skip this
      #weightmat <- weightmato/delta #as the delta changes, but not the weightmato
      #weightmat[!is.finite(weightmat)] <- 0
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
    weightmatorig <-stats::as.dist(weightmatorig) 
    #resmat <- weightmatm*as.matrix((deltam - doutm)^2)
    #spp <- colMeans(resmat) 
    weightmat <- stats::as.dist(weightmatm)
    spoint <- spp(delta, dout, weightmat)
    resmat<-spoint$resmat
    rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])
    spp <- spoint$spp
    if (verbose > 0 && itel == itmax) warning("Iteration limit reached! You may want to increase the itmax argument!")
    if (principal) {
        xnew_svd <- svd(xnew)
        xnew <- xnew %*% xnew_svd$v
    }
    #stressen <- sum(weightmat*(doute-delta)^2)
    if(verbose>1) cat("*** Stress:",snew, "; Stress 1 (default reported):",sqrt(snew),"\n")  
    out <- list(delta=deltaorig, dhat=delta, confdist=dout, iord=dhat2$iord.prim, conf = xnew, stress=sqrt(snew), spp=spp,  ndim=p, weightmat=weightmatorig, resmat=resmat, rss=rss, init=xstart, model="Sammon SMACOF", niter = itel, nobj = dim(xnew)[1], type = type, call=match.call(), stress.m=snew, alpha = anew, sigma = snew, tdelta=deltaold, parameters=c(kappa=1,lambda=1), pars=c(kappa=1,lambda=1), theta=c(kappa=1,lambda=1),tweightmat=weightmat)
    class(out) <- c("smacofP","smacofB","smacof")
    out
  }

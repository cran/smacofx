#' R stress SMACOF
#'
#' An implementation to minimize r-stress by majorization with ratio, interval, monotonic spline and ordinal optimal scaling. Uses a repeat loop.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param r power of the transformation of the fitted distances (corresponds to kappa/2 in power stress); defaults to 0.5 for standard stress
#' @param type what type of MDS to fit. Currently one of "ratio", "interval", "mspline" or "ordinal". Default is "ratio".
#' @param ties the handling of ties for ordinal (nonmetric) MDS. Possible are "primary" (default), "secondary" or "tertiary".
#' @param weightmat a matrix of finite weights. 
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should fitting information be printed; if > 0 then yes
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration
#' @param spline.degree Degree of the spline for ‘mspline’ MDS type
#' @param spline.intKnots Number of interior knots of the spline for ‘mspline’ MDS type
#'
#' @return a 'smacofP' object (inheriting from 'smacofB', see \code{\link[smacof]{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed, untransformed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Transformed fitted configuration distances
#' \item iord: Optimally scaled disparities function
#' \item conf: Matrix of fitted configuration
#' \item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
#' \item spp: Stress per point 
#' \item ndim: Number of dimensions
#' \item weightmat: Weighting matrix as supplied
#' \item resmat: Residual matrix
#' \item rss: Sum of residuals
#' \item init: The starting configuration
#' \item model: Name of MDS model
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item type: Type of optimal scaling 
#' \item call : the matched call
#' \item stress.m: Default stress (stress-1^2)
#' \item alpha: Alpha matrix
#' \item sigma: Stress
#' \item parameters, pars, theta: Optimal transformation parameter
#' \item tweightmat: Transformed weighting matrix (here NULL)
#'}
#'
#'
#' @importFrom stats dist as.dist
#' 
#' @seealso \code{\link[smacof]{smacofSym}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#'
#' ## ordinal MDS
#' res<-rStressMin(as.matrix(dis), type = "ordinal", r = 1, itmax = 1000)
#' res
#' summary(res)
#' plot(res)
#'
#' ## spline MDS 
#' ress<-rStressMin(as.matrix(dis), type = "mspline", r = 1,
#'       itmax = 1000)
#' ress
#' plot(ress,"Shepard")
#' 
#' @export
rStressMin <- function(delta, r=0.5, type=c("ratio","interval","ordinal","mspline"), ties="primary", weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal = FALSE, spline.degree = 2, spline.intKnots = 2) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    #if(missing(weightmat)) weightmat <-
    if(inherits(weightmat,"dist") || is.data.frame(weightmat)) weightmat <- as.matrix(weightmat)
    if(!isSymmetric(weightmat)) stop("weightmat is not symmetric.\n")
    #r <- kappa/2
    ## -- Setup for MDS type
    if(missing(type)) type <- "ratio"
    type <- match.arg(type, c("ratio", "interval", "ordinal","mspline"), several.ok = FALSE)
    trans <- type
    typo <- type
    if (trans=="ratio") {
        trans <- "none"
    }
    else if (trans=="ordinal" && ties=="primary") {
        trans <- "ordinalp"
        typo <- "ordinal (primary)"
        }
    else if(trans=="ordinal" && ties=="secondary") {
        trans <- "ordinals"
        typo <- "ordinal (secondary)"
        }
    else if(trans=="ordinal" && ties=="tertiary") {
        trans <- "ordinalt"
        typo <- "ordinal (tertiary)"
        }
    else if(trans=="spline"){
        trans <- "mspline"
        typo <- "mspline"
        type <- "mspline"
        }
    if(verbose>0) cat(paste("Minimizing",type,"rStress with r=",r,"\n"))
    n <- nrow (delta)
    normi <- 0.5
    ##normi <- n #if normi=n we can use the iord structure in plot.smacofP
    ## but the problem is we don't get the correct stress then anymore.
    p <- ndim
    if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
    if(is.null(rownames(delta))) rownames(delta) <- 1:n 
    labos <- rownames(delta) #labels
    deltaorig <- delta
    #delta <- delta^lambda
    #weightmat <- weightmat^nu
    weightmat[!is.finite(weightmat)] <- 0
    disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
    delta <- delta / enorm (delta, weightmat) #CHECK: Was before before disobj before
    ## Add an intercept to the spline base transformation
    if (trans == "mspline") disobj$base <- cbind(rep(1, nrow(disobj$base)), disobj$base)
    deltaold <- delta
    itel <- 1
    ##Starting Configs
    xold  <- init
    if(is.null(init)) xold <- smacof::torgerson (delta, p = p)
    xstart <- xold
    xold <- xold / enorm(xold) 
    nn <- diag(n)
    dold <- sqdist(xold)
   ##first optimal scaling
    ## eold <- as.dist(sqrt(dold)) #was bug prior to 1.6-1
    eold <- as.dist(mkPower(dold,r))
    dhat <- smacof::transform(eold, disobj, w = as.dist(weightmat), normq = normi)
    dhatt <- dhat$res 
    dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
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
      #e <- as.dist(sqrt(dnew)) #I need the dist(x) here for interval, #was bug prior to 1.6-1
      e <- as.dist(mkPower(dnew,r)) 
      dhat2 <- smacof::transform(e, disobj, w = as.dist(weightmat), normq = normi)  ## dhat update
      dhatt <- dhat2$res 
      dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
      delta <- as.matrix(dhatd)
      rnew <- sum (weightmat * delta * mkPower (dnew, r))
      nnew <- sum (weightmat * mkPower (dnew, 2 * r))
      anew <- rnew / nnew
      snew <- 1 - 2 * anew * rnew + (anew ^ 2) * nnew
      if(is.na(snew)) { #if there are issues in the values
          snew <- sold
          dnew <- dold
          anew <- aold
          xnew <- xold
      }   
      if (verbose>2) {
        cat(
          formatC(itel, width = 4, format = "d"),
          formatC(
            sold,
            digits = 10,
            width = 13,
            format = "f"
          ),
          formatC(
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
    if (verbose > 0 && itel == itmax) warning("Iteration limit reached! You may want to increase the itmax argument!")
    if (principal) {
        xnew_svd <- svd(xnew)
        xnew <- xnew %*% xnew_svd$v
    }
    #stressen <- sum(weightmat*(doute-delta)^2)
    if(verbose>1) cat("*** Stress:",snew, "; Stress 1 (default reported):",sqrt(snew),"\n")
    #delta is input delta, tdelta is input delta with explicit transformation and normalized, dhat is dhats 
    out <- list(delta=deltaorig, tdelta=deltaold, dhat=delta, confdist=dout, iord=dhat2$iord.prim, conf = xnew, stress=sqrt(snew), spp=spp,  ndim=p, weightmat=weightmat, resmat=resmat, rss=rss, init=xstart, model="r-stress MDS", niter = itel, nobj = dim(xnew)[1], type = type, call=match.call(), stress.m = snew, alpha = anew, sigma = snew, parameters=c(r=r), pars=c(r=r), theta=c(r=r), tweightmat=NULL)
    class(out) <- c("smacofP","smacofB","smacof")
    out
  }


#' @rdname rStressMin
#' @export
rstressMin <- rStressMin

#' @rdname rStressMin
#' @export
rstressmds <- rStressMin

#' @rdname rStressMin
#' @export
rstress <- rStressMin




#' Power Stress SMACOF
#'
#' An implementation to minimize power stress by majorization with ratio or interval optimal scaling. Usually more accurate but slower than powerStressFast. Uses a repeat loop.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param kappa power of the transformation of the fitted distances; defaults to 1
#' @param lambda the power of the transformation of the proximities; defaults to 1
#' @param nu the power of the transformation for weightmat; defaults to 1
#' @param type what type of MDS to fit. One of "ratio" or "interval". Default is "ratio".
#' @param weightmat a matrix of finite weights or dist object
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; if > 1 then yes
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration
#'
#' @return a 'smacofP' object (inheriting from 'smacofB', see \code{\link[smacof]{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed, untransformed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Transformed fitted configuration distances
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
#' \item tweightmat: transformed weighthingmatrix (here weightmat^nu)
#' }
#'
#' @importFrom stats dist as.dist
#' 
#' @seealso \code{\link[smacof]{smacofSym}}
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-powerStressMin(dis,type="ratio",kappa=2,lambda=1.5,itmax=1000)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
powerStressMin <- function (delta, kappa=1, lambda=1, nu=1, type="ratio", weightmat=1-diag(nrow(delta)), init=NULL, ndim = 2, acc= 1e-6, itmax = 10000, verbose = FALSE, principal=FALSE) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
    #if(missing(weightmat))  weightmat <- 
    if(inherits(weightmat,"dist") || is.data.frame(weightmat)) weightmat <- as.matrix(weightmat)
    if(!isSymmetric(weightmat)) stop("weightmat is not symmetric.\n")
    type <- match.arg(type, c("ratio", "interval"),several.ok = FALSE) 
    trans <- type
    typo <- type
    if (trans=="ratio"){
    trans <- "none"
    }
    if(verbose>0) cat("Minimizing",type,"power-stress with kappa=",kappa,"lambda=",lambda,"nu=",nu,"\n")

    n <- nrow (delta)
    r <- kappa/2
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
    disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)
    deltaold <- delta
    itel <- 1
    xold <- init
    if(is.null(init)) xold <- smacof::torgerson (delta, p = p)
    xstart <- xold
    xold <- xold / enorm (xold)
    n <- nrow (xold)
    nn <- diag (n)
    dold <- sqdist (xold)
    #Optimal scaling
    #eold <- as.dist(sqrt(dold)) #was bug prior to 1.6-1
    eold <- as.dist(mkPower(dold,r))
    dhat <- smacof::transform(eold, disobj, w = as.dist(weightmat), normq = 0.5)
    dhatt <- dhat$res 
    dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
    delta <- as.matrix(dhatd)
    rold <- sum (weightmat * delta * mkPower (dold, r))
    nold <- sum (weightmat * mkPower (dold, 2 * r))
    aold <- rold / nold
    sold <- 1 - 2 * aold * rold + (aold ^ 2) * nold
    #Optimizing
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
      #Optimal scaling 
      #e <- as.dist(sqrt(dnew)) #was bug prior to 1.6-1
      e <- as.dist(mkPower(dnew, r))
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
     }
     xnew <- xnew/enorm(xnew)
     ## relabeling as they were removed in the optimal scaling
     rownames(delta) <- labos
     attr(xnew,"dimnames")[[1]] <- rownames(delta)
     attr(xnew,"dimnames")[[2]] <- paste("D",1:p,sep="")
     doutm <- mkPower(sqdist(xnew),r)
     deltam <- delta
     delta <- stats::as.dist(delta)
     #delta <- structure(delta, Size = n, call = quote(as.dist.default(m=b)),
     #                  class = "dist", Diag = FALSE, Upper = FALSE)
     deltaorig <- stats::as.dist(deltaorig)
     deltaold <- stats::as.dist(deltaold)
     #doute <- doutm/enorm(doutm) #this is an issue here!
     #doute <- stats::as.dist(doute)
     dout <- stats::as.dist(doutm)
     weightmatm <-weightmat
     weightmat <- stats::as.dist(weightmatm)
     weightmato <- stats::as.dist(weightmato)
     #resmat <- weightmatm*as.matrix((deltam - doutm)^2) #old spp now we sum up to 100
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
     if(verbose>1) cat("*** Stress:",snew, "; Stress 1 (default reported):",sqrt(snew), "\n")  
    out <- list(delta=deltaorig, dhat=delta, confdist=dout, iord=dhat2$iord.prim, conf = xnew, stress=sqrt(snew), spp=spp,  ndim=p, weightmat=weightmato, resmat=resmat, rss=rss, init=xstart, model="Power-Stress SMACOF", niter = itel,nobj = dim(xnew)[1], type = type, call=match.call(), stress.m=snew, alpha = anew, sigma = snew, tdelta=deltaold, parameters=c(kappa=kappa,lambda=lambda,nu=nu), pars=c(kappa=kappa,lambda=lambda,nu=nu), theta=c(kappa=kappa,lambda=lambda,nu=nu),tweightmat=weightmat)
    class(out) <- c("smacofP","smacofB","smacof")
    out
}


#' @rdname powerStressMin
#' @export
powerstressMin <- powerStressMin

#' @rdname powerStressMin
#' @export
postmds <- powerStressMin

#' @rdname powerStressMin
#' @export
pstressMin <- powerStressMin

#' @rdname powerStressMin
#' @export
pStressMin <- powerStressMin

#' @rdname powerStressMin
#' @export
pstressmds <- powerStressMin

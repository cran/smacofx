#' Nonlinear ratio MDS with optimal power of dissimilarities
#'
#' An implementation to minimize explicitly normalized stress over dissimilarities to a power by majorization with ratio optimal scaling in an alternating minimization algorithm. The optimal power transformation lambda of the dissimilarities is found by an inner optimization step via the Brent-Dekker method.
#' 
#' @param delta dist object or a symmetric, numeric data.frame or matrix of distances
#' @param type what type of MDS to fit. Currently only "ratio".
#' @param weightmat a matrix of finite weights. 
#' @param init starting configuration
#' @param ndim dimension of the configuration; defaults to 2
#' @param acc numeric accuracy of the iteration. Default is 1e-6.
#' @param itmax maximum number of iterations. Default is 10000.
#' @param verbose should iteration output be printed; defaults to 'FALSE'.
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration.
#' @param interval the line constraints c(upper, lower), within which to look for the optimal power transformation lambda. Defaults to c(0,4).
#'
#' @return a 'smacofP' object (inheriting from 'smacofB', see \code{\link[smacof]{smacofSym}}). It is a list with the components
#' \itemize{
#' \item delta: Observed, untransformed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
#' \item confdist: Transformed fitted configuration distances
#' \item iord: optimal scaling ordering
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
#' \item tweightmat: transformed weighting matrix (here NULL)
#' \item pars, theta: The optimal transformation parameter lambda
#'}
#'
#'
#' @importFrom stats dist as.dist optimize
#' 
#' @seealso See \code{\link[stops]{stops}} for a similar, more flexible idea.
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-opmds(dis,itmax=1000)
#' res
#' summary(res)
#' plot(res)
#'
#' @export
opmds <- function(delta, type="ratio", weightmat=1-diag(nrow(delta)), init=NULL, ndim=2, itmax = 1000, acc = 1e-10, verbose = FALSE, principal=FALSE, interval = c(0, 4)) {
    if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
    if(!isSymmetric(delta)) stop("delta is not symmetric.\n")
    if(inherits(weightmat,"dist") || is.data.frame(weightmat)) weightmat <- as.matrix(weightmat)
    if(!isSymmetric(weightmat)) stop("weightmat is not symmetric.\n")
    ## -- Setup for MDS type
    if(missing(type)) type <- "ratio"
    type <- match.arg(type, c("ratio"))
    trans <- "none"
    nobj <- nrow(delta)
    if (ndim > (nobj - 1)) stop("Maximum number of dimensions is n-1!")
    if(is.null(rownames(delta))) rownames(delta) <- 1:nobj
    labos <- rownames(delta) #labels
    deltaorig <- delta
    weightmat[!is.finite(weightmat)] <- 0
    if (interval[1] == interval[2]) {
      r <- interval[1]
      fixed <- TRUE
    } else {
      r <- (interval[1] + interval[2]) / 2
      fixed <- FALSE
    }
    ep <- delta ^ r
    ep <- ep / enorm (ep, weightmat)
    #disobj <- smacof::transPrep(as.dist(delta), trans = trans, spline.intKnots = 2, spline.degree = 2)
    deltaold <- ep
    xe <- init
    if (is.null(xe)) {
      dd <- ep ^ 2 
      rd <- rowSums(dd) / nobj
      sd <- mean(ep)
      ce <- -0.5 * (dd - outer(rd, rd) + sd)
      ee <- eigen(ce)
      v <- pmax(ee$values, 0)
      if (ndim == 1) 
           normdiag <- cbind(sqrt(v[1]))
      else normdiag <- diag(sqrt(v[1:ndim]))
      xe <- ee$vectors[, 1:ndim] %*% normdiag
    }
    if(ncol(xe) != ndim) stop("Column number of init matrix and ndim argument must be the same.")
    rownames(xe) <- labos
    de <- as.matrix(dist(xe))
    g <- function(r, delta, de, weightmat) {
        ##adapted this to enorm the delta
        ep <- delta ^ r
        ep <- ep / enorm(ep, weightmat)
        return(sum(weightmat*(ep - de) ^ 2))
    }
    sold <- sum(weightmat*(ep - de) ^ 2)
    itel <- 1
    repeat {
      b <- -ep * ifelse(de == 0, 0, 1 / de)
      diag(b) <- -rowSums(b)
      xe <- (b %*% xe) / nobj
      de <- as.matrix(dist(xe))
      ##not sure if we need this optimal scaling as we just do ratio: iord should then simply be 1:nobj 
#     dhat2 <- smacof::transform(dist(xe), disobj, w = as.dist(weightmat), normq = 0.5)  ## dhat update
      if (!fixed) {
          r <- stats::optimize(g, interval = interval, delta = deltaorig, de = de, weightmat = weightmat)$minimum
          ##the interval argument are the box constrains for optimize(), delta, de and weightmat are extra args for g()
      }
      ep <- delta ^ r
      ep <- ep / enorm(ep, weightmat)
#      dhatt <- dhat2$res 
#      dhatd <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
#     ep <- as.matrix(dhatd)
      snew <- sum(weightmat*(ep - de) ^ 2)
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(sold, digits = 6, format = "f"),
          "snew ",
          formatC(snew, digits = 6, format = "f"),
          "pow  ",
          formatC(r, digits = 6, format = "f"),
          "\n"
        )
      }
      if (((sold - snew) < acc) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      sold <- snew
    }
    if (principal) {
        xe_svd <- svd(xe)
        xe <- xe %*% xe_svd$v
    }
      out <- list(delta=deltaorig,
      tdelta = deltaold,
      dhat=ep, 
      confdist = de,
      iord=order(as.vector(delta)),
      conf = xe,
      stress = sqrt(snew),
      spp=NULL,
      ndim=ndim,
      model= c("nonlinear MDS with Power"),
      niter = itel,
      nobj=nobj,
      type=type,
      weightmat=weightmat,
      stress.m=snew,
      tweightmat=NULL,
      parameters=c(lambda=r),
      pars=c(lambda=r),
      theta=c(lambda=r),
      call=match.call())
      class(out)<-c("smacofP","smacofB")
      return(out)
}

## smacofPO <-
##   function(delta,
##            interval = c(0, 4),
##            xe = NULL,
##            itmax = 1000,
##            eps = 1e-10,
##            verbose = FALSE) {
##     nobj <- nrow(delta)
##     if (is.null(xe)) {
##       dd <- delta ^ 2
##       rd <- rowSums(dd) / nobj
##       sd <- mean(delta)
##       ce <- -.5 * (dd - outer(rd, rd) + sd)
##       ee <- eigen(ce)
##       xe <- ee$vectors[, 1:2] %*% diag(sqrt(ee$values[1:2]))
##     }
##     de <- as.matrix(dist(xe))
##     if (interval[1] == interval[2]) {
##       r <- interval[1]
##       fixed <- TRUE
##     } else {
##       r <- (interval[1] + interval[2]) / 2
##       fixed <- FALSE
##     }
##     g <- function(r, delta, de) {
##       return(sum(((delta ^ r) - de) ^ 2))
##     }
##     ep <- delta ^ r
##     sold <- sum((ep - de) ^ 2)
##     itel <- 1
##     repeat {
##       b <- -ep * ifelse(de == 0, 0, 1 / de)
##       diag(b) <- -rowSums(b)
##       xe <- (b %*% xe) / nobj
##       de <- as.matrix(dist(xe))
##       smid <- sum((ep - de) ^ 2)
##       if (!fixed) {
##         r <- optimize(g, interval = interval, delta = delta, de = de)$minimum
##       }
##       ep <- delta ^ r
##       snew <- sum((ep - de) ^ 2)
##       if (verbose) {
##         cat(
##           "itel ",
##           formatC(itel, format = "d"),
##           "sold ",
##           formatC(sold, digits = 6, format = "f"),
##           "smid ",
##           formatC(smid, digits = 6, format = "f"),
##           "snew ",
##           formatC(snew, digits = 6, format = "f"),
##           "pow  ",
##           formatC(r, digits = 6, format = "f"),
##           "\n"
##         )
##       }
##       if (((sold - snew) < 1e-10) || (itel == itmax)) {
##         break
##       }
##       itel <- itel + 1
##       sold <- snew
##     }
##     return(list(
##       x = xe,
##       d = de,
##       e = ep,
##       r = r,
##       itel = itel,
##       stress.m = snew,
##       stress= sqrt(snew/sum(ep^2)) 
##     ))
##   }

#' MDS Jackknife for smacofP objects
#'
#' These functions perform an MDS Jackknife and plot the corresponding solution.
#'
#' @param object  Object of class smacofP if used as method or another object inheriting from smacofB (needs to be called directly as jackmds.smacofP then).
#' @param eps Convergence criterion
#' @param itmax Maximum number of iterations 
#' @param verbose If 'TRUE', intermediate stress is printed out.
#'
#' @details  In order to examine the stability solution of an MDS, a Jackknife on the configurations can be performed (see de Leeuw & Meulman, 1986) and plotted. The plot shows the jackknife configurations which are connected to their centroid. In addition, the original configuration (transformed through Procrustes) is plotted. The Jackknife function itself returns also a stability measure (as ratio of between and total variance), a measure for cross validity, and the dispersion around the original smacof solution.
#'
#' @return An object of class 'smacofJK', see \code{\link[smacof]{jackmds}}. With values 
#' \itemize{
#' \item smacof.conf: Original configuration
#' \item jackknife.confboot: An array of n-1 configuration matrices for each Jackknife MDS solution
#' \item comparison.conf: Centroid Jackknife configurations (comparison matrix)
#' \item cross: Cross validity
#' \item stab: Stability coefficient
#' \item disp: Dispersion
#' \item loss: Value of the loss function (just used internally)
#' \item ndim: Number of dimensions
#' \item call: Model call
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' }
#'
#'
#' @importFrom smacof jackmds
#' 
#' @export
#' @examples
#' dats <- na.omit(smacof::PVQ40[,1:5])
#' diss <- dist(t(dats))   ## Euclidean distances
#' fit <- rStressMin(diss,type="ordinal",r=0.4,itmax=1000) ## 2D ordinal MDS
#'
#' res.jk <- jackmds(fit)
#' plot(res.jk, col.p = "black", col.l = "gray")
#' plot(res.jk, hclpar = list(c = 80, l = 40))
#' plot(res.jk, hclpar = list(c = 80, l = 40), plot.lines = FALSE)

jackmds.smacofP <- function(object, eps = 1e-6, itmax = 100, verbose = FALSE) 
{
  ## object... object of class smacofP
  #if (class(object)[1] != "smacofP") stop("Jackknife is currently implemented for objects of class smacofB from smacofSym() only! \n")
  ##  if (object$model == "SMACOF constraint") stop("Jackknife is currently implemented for smacofSym() objects only! \n")

    calli <- match.call()
    delta <- as.matrix(object$delta)
    n <- nrow(delta)
    x0 <- object$conf
    ndim <- object$ndim
    type <- object$type
    weightmat <- as.matrix(object$weightmat)
    init <- as.matrix(object$init)
    if(!is.numeric(object$call$itmax)) itmaxi <- object$niter+1000 else itmaxi <- object$call$itmax #this is for using it in the jackmds.stops method. We have not substituted itmaxi there and it would be hard to do, so I just use the returned iterations upon convergence and add 1000 iterations as max in case it was passed down. Shouldn't effect a direct smacofP object. 
  
    
   #x0 <- smacofSym (delta, ndim = ndim, metric = metric) $ conf
    tmpo <- smacofxDeleteOne(object, delta, weightmat = weightmat, init=init, ndim = ndim, type = type, verbose=verbose,itmaxi=itmaxi)
    xx <- tmpo$x
    kk <- array(rep(diag (ndim), n), c(ndim, ndim, n))
    cc <- matrix(0, n, ndim)
    bb <- matrix(0, n, ndim)
    yy <- xx
    oloss <- Inf
    itel <- 1
    repeat {
        y0 <- matrix (0, n, ndim)
        for (i in 1:n) {
            y0 <- y0 + xx[, , i] %*% kk[, , i]
            }
        y0 <- ((n - 1) * y0) / (n * (n - 2))
        for (i in 1:n) {
            zz <- matrix (0, n, ndim)
            for (j in i : n) {
                zz <- zz + xx[, , j] %*% kk[, , j]
                }
            xz <- crossprod ( xx[, , i], zz)
            kk[, , i] <- procruster (xz)
            }
        nloss <- 0
        for (i in 1:n) {
            yy[, , i] <- xx[, , i] %*% kk [, , i] 
            yy[i, , i] <- n * y0[i, ] / (n - 1)
            yy[, , i] <- yy[, , i] - outer (rep (1, n), y0[i, ] / (n - 1))
            nloss <- nloss + sum ( (y0 - yy[, , i]) ^ 2)
            }
        if (verbose) {
  	     cat("Iteration: ", formatC (itel, digits=3, width=3), "Old Loss: ", formatC (oloss, digits=10, width=15, format="f"), 
  	         "New Loss: ", formatC (nloss, digits=10, width=15, format="f"), "\n")
			 }
	    if (((oloss - nloss) < eps) || (itel == itmax)) {
	        break()
	        }
	    itel <- itel + 1
	    oloss <- nloss
	  }
	 x0 <- x0 %*% procruster(crossprod (x0, y0))
   
   if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")
	
  ## stability measure
  stab.num <- sum(apply(yy, 3, function(yyi) (norm(yyi-y0))^2))  ## numerator stab
  stab.denom <- sum(apply(yy, 3, function(yyi) (norm(yyi))^2))   ## denominator stab
  stab <- 1 - stab.num/stab.denom

  ## cross vallidity
  cross.num <- n*(norm(x0-y0))^2
  cross <- 1 - cross.num/stab.denom
  
  ## dispersion around x0
  disp <- 2 - (stab + cross)

  result <- list(smacof.conf = x0, jackknife.conf = yy, comparison.conf = y0, stab = stab, cross = cross, disp = disp, niter = itel, loss = nloss, nobj = n, call = calli)
  class(result) <- "smacofJK"
  result
}


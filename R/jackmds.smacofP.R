#' MDS Jackknife for smacofP objects
#'
#' These methods perform an MDSJackknife and plot the corresponding solution.
#'
#' @param object  Object of class smacofP if used as method or another object inheriting from smacofB (needs to be called directly as bootmds.smacofP then).
#' @param eps Convergence criterion
#' @param itmax Maximum number of iterations 
#' @param verbose If 'TRUE', intermediate stress is printed out.
#'
#' @details  In order to examine the stability solution of an MDS, a Jackknife on the configurations can be performed (see de Leeuw & Meulman, 1986) and plotted. The plot shows the jackknife configurations which are connected to their centroid. In addition, the original configuration (transformed through Procrustes) is plotted. The Jackknife function itself returns also a stability measure (as ratio of between and total variance), a measure for cross validity, and the dispersion around the original smacof solution.
#'
#' @return An object of class 'smacofJK', see \code{\link[smacof]{jackmds}}. With values 
#' \itemize{
#' \item smacof.conf: Configurations
#' \item jackknife.confboot: An array of n-1 configuration matrices for each Jackknife MDS solution
#' \item comparison.conf: Centroid Jackknife configurations (comparison matrix)
#' \item cross: Cross validity
#' \item stab: Stability coefficient
#' \item disp: Dispersion
#' \item loss: Value of the loss function
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
#' data <- na.omit(smacof::PVQ40[,1:5])
#' diss <- dist(t(data))   ## Euclidean distances
#' fit <- rStressMin(diss,"ratio",r=0.5,itmax=1000) ## 2D ratio MDS
#'
#' res.jk <- jackmds(fit)
#' plot(res.jk, col.p = "black", col.l = "gray")
#' plot(res.jk, hclpar = list(c = 80, l = 40))
#' plot(res.jk, hclpar = list(c = 80, l = 40), plot.lines = FALSE)
#' 
jackmds.smacofP <- function(object, eps=1e-06,itmax=100,verbose = FALSE)
    #TODO add an itmax argument?
{
    if(any(class(object)=="smacofB"))
    {
        class(object) <- c("smacofB",class(object))
    } else stop("MDS object must inherit from smacofB.")
    out <- smacof::jackmds(object,eps=eps,itmax=itmax,verbose=verbose)
    class(object) <- class(object)[-1]
    out
 }

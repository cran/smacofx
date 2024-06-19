#' Exploring initial configurations in an agnostic way  
#'
#' Allows to user to explore the effect of various starting configurations when fitting an MDS model. This is a bit more general than the icExplore function in smacof, as we allow any PS model to be used as the model is either setup by call or by a prefitted object (for the models in cops and stops we do not have a single UI function which necessitates this). Additionally, one can supply their own configurations and not just random ones.    
#'
#' @param object A fitted object of class 'smacofP', 'smacofB' or 'smacof'. If supplied this takes precedence over the call argument. If given this is added to the output and may be the optimal one. 
#' @param mdscall Alternatively to a fitted object, one can pass a syntactically valid call for any of the MDS functions cops, stops or smacof that find a configuration (not the ones that do parameter selection like pcops or stops). If object and call is given, object takes precedence.
#' @param ndim Dimensions of target space.
#' @param conflist Optional list of starting configurations.
#' @param nrep If conflist is not supplied, how many random starting configurations should be used. 
#' @param returnfit Should all fitted MDS be returned. If FALSE (default) none is returned.
#' @param verbose If >0 prints the fitting progress.
#' @param min lower bound for the uniform distribution to sample from
#' @param max upper bound for the uniform distribution to sample from
#'
#' @details  
#' If no configuration list is supplied, then nrep configurations are simulated. They are drawn from a ndim-dimensional uniform distribution with minimum min and maximum max. We recommend to use the route via supplying a fitted model as these are typically starting from a Torgerson configuration and are likely quite good.    
#'
#' @export
#' @return an object of class 'icexplore', see \code{\link[smacof]{icExplore}} for more. There is a plot method in package 'smacof'. 
#'
#' @importFrom smacof Procrustes sim2diss mds
#' @importFrom stats cor runif
#' 
#' @examples
#'
#' dis<-kinshipdelta
#' 
#' ## Version 1: Using a fitted object (recommended)
#' res1<-rStressMin(dis,type="ordinal",itmax=100)
#' resm<-icExploreGen(res1,nrep=5)
#'
#' ## Version 2: Using a call object and supplying conflist
#' conflist<-list(res1$init,jitter(res1$init,1),jitter(res1$init,1),jitter(res1$init,1))
#' c1 <- call("smds",delta=dis,tau=0.2,itmax=100)
#' resm<-icExploreGen(mdscall=c1,conflist=conflist,returnfit=TRUE)
#'
#' plot(resm)
#' 
icExploreGen <- function(object, mdscall=NULL, conflist, nrep = 100, ndim, returnfit = FALSE, min=-5, max=5, verbose=FALSE) 
{
    #diss <- delta
    #if ((is.matrix(diss)) || (is.data.frame(diss))) {
    #    diss <- strucprep(diss)
    #    attr(diss, "Labels") <- rownames(delta)
    #}
    #checkdiss(diss)
    #n <- attr(diss, "Size")
    if(missing(ndim)) ndim <- 2
    if(!missing(conflist)) nrep <- length(conflist)
    if(is.null(mdscall) && missing(object)) stop("An MDS function call (either as call object or list) or a fitted object must be provided.")
    if(is.list(mdscall)) mdscall <- as.call(mdscall)
    if(!missing(object))
    {
        n <-object$nobj
        mdscall <- object$call
        ndim <- object$ndim
    } else n <- nrow(as.matrix(mdscall$delta))
    v.stress <- vector()
    configs <- list()
    if(missing(conflist)) conflist <- replicate(nrep, matrix(stats::runif(n * ndim, min = min, max = max), nrow = n, ncol = ndim),simplify=FALSE) 
    simi <- matrix(0, nrow = nrep, ncol = nrep)
    labels <- as.character(1:nrep)
    restot <- list()
    for (i in 1:nrep) {
        if (verbose) 
            cat("IC:", i, "\n")
        mdscall$init <- conflist[[i]]
        aus1 <- eval(mdscall)
        configs[[i]] <- aus1$conf
        v.stress[i] <- aus1$stress
        if(length(grep("cop",as.character(mdscall[[1]])))>0) v.stress[i] <- aus1$copstress #check for cops
        restot[[i]] <- aus1
    }
    for (i in 2:nrep) {
        je <- i - 1
        for (j in 1:je) {
            a <- configs[[i]]
            b <- configs[[j]]
            aus <- smacof::Procrustes(a, b)
            a1 <- c(a)
            b1 <- c(aus$Yhat)
            simi[i, j] <- stats::cor(a1, b1)
        }
    }
    if(is.null(mdscall$itmax)) itmax <- 10000 else itmax <- mdscall$itmax
    sim2 <- simi + t(simi) + diag(nrep)
    diss1 <- smacof::sim2diss(sim2, method = "corr")
    aus2 <- smacof::mds(diss1, type = "interval", itmax = itmax)
    x <- aus2$conf[, 1]
    y <- aus2$conf[, 2]
    if (!returnfit) 
        restot <- NULL
    res <- list(mdsfit = restot, stressvec = v.stress, conf = cbind(x, 
        y), nrep = nrep, call = match.call())
    class(res) <- "icexplore"
    return(res)
}

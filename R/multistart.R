#' Multistart MDS function
#'
#' For different starting configurations, this function fits a series of PS models given in object or call and returns the one with the lowest stress overall. The starting configuirations can be supplied or are generated internally.
#'
#' @param object A fitted object of class 'smacofP', 'smacofB' or 'smacof'. If supplied this takes precedence over the call argument. If given this is added to the output and may be the optimal one. 
#' @param mdscall Alternatively to a fitted object, one can pass a syntactically valid call for any of the MDS functions cops, stops or smacof that find a configuration (not the ones that do parameter selection like pcops or stops). If object and call is given, object takes precedence.
#' @param ndim Dimensions of target space.
#' @param conflist Optional list of starting configurations.
#' @param nstarts If conflist is not supplied, how many random starting configurations should be used. The default is 108, which implies that at least one of the stress is within the lowest 1 percent of all stresses with probability of 1/3 or within the lowest 5 percent of stresses with probability 0.996  
#' @param return.all Should all fitted MDS be returned. If FALSE (default) only the optimal one is returned.
#' @param verbose If >0 prints the fitting progress.
#' @param min lower bound for the uniform distribution to sample from
#' @param max upper bound for the uniform distribution to sample from
#'
#' @details  
#' If no configuration list is supplied, then nstarts configurations are simulated. They are drawn from a ndim-dimesnional uniform distribution with minimum min and maximum max. We recommend to use the route via supplying a fitted model as these are typically starting from a Torgerson configuration and are likely quite good.    
#'
#' One can simply extract $best and save that and work with it right away.
#' 
#' @return if 'return.all=FALSE', a list with the best fitted model as '$best' (minimal badness-of-fit of all fitted models) and '$stressvec' the stresses of all models. If 'return.all=TRUE' a list with slots
#' \itemize{
#' \item best: The object resulting from the fit that had the overall lowest objective function value (usually stress)
#' \item stressvec: The vector of objective function values
#' \item models: A list of all the fitted objects.
#' }
#'
#' @importFrom stats runif
#' 
#' @export
#' @examples
#' dis<-smacof::kinshipdelta
#' 
#' ## Version 1: Using a fitted object (recommended)
#' res1<-rStressMin(delta=dis,type="ordinal",itmax=100)
#' resm<-multistart(res1,nstarts=2)
#' ## best model
#' res2<-resm$best
#' #it's starting configuration
#' res2$init
#'
#' ## Version 2: Using a call object and supplying conflist
#' conflist<-list(res2$init,jitter(res2$init,1))
#' c1 <- call("rstressMin",delta=dis,type="ordinal",itmax=100)
#' resm<-multistart(mdscall=c1,conflist=conflist,return.all=TRUE)
#' 
multistart <- function(object,mdscall=NULL,ndim=2,conflist,nstarts=108,return.all=FALSE,verbose=TRUE,min=-5, max=5)
{
    # We need the call because there isn't a singular or nested API for all the functions and it would be undoable/annoying to build all the permutations into the function
    if(is.null(mdscall) && missing(object)) stop("An MDS function call (either as call object or list) or a fitted object must be provided.") 
    if(is.call(mdscall))
    {
        delta <- mdscall$delta
        n <- nrow(as.matrix(delta))
    }
    if(is.list(mdscall))
    {
        mdscall <- as.call(mdscall)
        delta <- mdscall$delta
        n <- nrow(as.matrix(delta))
    } 
    if(!missing(object))
        {
            mdscall <- object$call
            delta <- object$delta
            n <- object$nobj
        }
    if(missing(conflist))
    {
        conflist <- replicate(nstarts,matrix(stats::runif(n*ndim,min=min,max=max),ncol=ndim),simplify=FALSE)
    }
    o <- vector("list",length(conflist))
    for(i in 1:length(conflist))
    {
         if(verbose>0) cat("Fitting Model:", i)
         mdscall$init <- conflist[[i]]
         o[[i]] <- eval(mdscall)
         if(verbose>1) cat(" Stress: ", o[[i]]$stress)
         if(verbose>0) cat("\n")
    }
    if(!missing(object)) o <- c(o,list(object)) #add the supplied object to the list
    stressvec <- sapply(o,function(x) x$stress)
    #if(!missing(object)) #we need to handle copsc differently as that one has copstress as the objective value, not stress.  
    #    {
    #        if(any(class(object)=="copsc"))  stressvec <- sapply(o,function(x) x$copstress)
    #    }
    if(length(grep("cop",as.character(mdscall[[1]])))>0) stressvec <- sapply(o,function(x) x$copstress)
    optimal <- which.min(stressvec)
    #TODO: if we want to return only an object
    #if(!return.all)
    #{
    #    return(o[[optimal]])
    #    } else {
    out <- list()
    out$best<- o[[optimal]]
    out$stressvec <- stressvec
    if(return.all) out <- c(out,models=list(o))
    return(out)
  }


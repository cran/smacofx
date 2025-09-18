#' Function to plot Procrustes aligned MDS configurations 
#'
#' This plots a series of Procrustes aligned configuration plots, aligned to the reference configuration. 
#' 
#'@param objectlist a list of object of class smacofP or smacofB (or anything with a $conf slot where the configuration is). Thse are the plots to Procrustes adjust.
#'@param reference the reference configuration/matrix
#'@param layoutmat a matrix layout for the plot region. It is the 'mat' argument to 'layout'. If missing, the function attempts to create a layout that can hold all plots, but it is not very smart. 
#'@param mains a character string of headings for each plot. It must be of length(objectlist).  
#'@param ... additional arguments passed to the confplot of 'plot.smacofP'. See 'plot.smacof' for the options.
#'
#' 
#'@importFrom graphics layout
#'@importFrom smacof Procrustes
#' 
#'@return no return value; just plots
#' 
#'@export
#' 
#'@examples
#'dis<-as.matrix(smacof::kinshipdelta)
#'res1<-powerStressMin(dis)
#'res0<-smacof::mds(dis)
#'res2<-smacof::mds(dis,type="ordinal")
#'alignplot(list(res1,res2),reference=res0$conf,mains=c("pstress","ordinal"),layoutmat=c(1,2))
alignplot <- function(objectlist, reference, layoutmat=NULL,mains=NULL,...)
{
   op <- par(no.readonly=TRUE)
   # if(length(objectlist)==1)
   # {
   #     if(is.null(mains)) mains <- "Configuration Plot"
   #     plot(objectlist[[1]],main=mains[1],...)
   #     
   # } else {
    if(is.null(mains)) mains <- rep("Configuration Plot",length(objectlist))
    if(!isTRUE(all.equal(length(mains),length(objectlist)))) stop("Argument mains must be of the same length as objectlist!")
    if(is.null(layoutmat))
    {
        plotnrs <- length(objectlist)
        if(length(objectlist)==1) matlt <- 1 else matlt <- matrix(seq(1,3*ceiling(plotnrs/3)),ncol=3,nrow=ceiling(plotnrs/3))
        layoutmat <- matlt
        }
        graphics::layout(mat=layoutmat)
        for(i in seq(1,length(objectlist))){
            objectlist[[i]]$conf <- smacof::Procrustes(reference,objectlist[[i]]$conf)$Yhat
            plot(objectlist[[i]],main=mains[[i]],...)
        }
    on.exit(par(op))
}

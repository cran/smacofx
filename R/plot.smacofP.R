#'S3 plot method for smacofP objects
#' 
#'@param x an object of class smacofP 
#'@param plot.type String indicating which type of plot to be produced: "confplot", "resplot", "Shepard", "stressplot","transplot", "bubbleplot" (see details)
#'@param plot.dim  dimensions to be plotted in confplot; defaults to c(1, 2)
#'@param main plot title
#'@param xlab label of x axis
#'@param ylab label of y axis
#'@param xlim scale of x axis
#'@param ylim scale of y axis
#'@param col vector of colors for the points
#'@param bubscale Scaling factor (size) for the bubble plot
#'@param label.conf List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color)
#'@param identify If 'TRUE', the 'identify()' function is called internally that allows to add configuration labels by mouse click
#'@param type What type of plot should be drawn (see also 'plot')
#'@param legend Flag whether legends should be drawn for plots that have legends
#'@param legpos Position of legend in plots with legends 
#'@param pch  Plot symbol
#'@param asp  Aspect ratio; defaults to 1 so distances between x and y are represented accurately; can lead to slighlty weird looking plots if the variance on one axis is much smaller than on the other axis; use NA if the standard type of R plot is wanted where the ylim and xlim arguments define the aspect ratio - but then the distances seen are no longer accurate
#'@param loess if TRUE a loess fit (by Tukey's rescending M-Estimator) of configuration distances explained by delta  is added to the Shepard plot
#'@param hull.conf Option to add convex hulls to a configuration plot. Hull index needs to be provided.
#'@param shepard.x Shepard plot only: original data (e.g. correlation matrix) can be provided for plotting on x-axis
#'@param cex Symbol size.
#'@param col.hist Color of the borders of the histogram.
#'@param shepard.lin Shepard plot only: if TRUE the Shepard plot is linearized so d^kappa~delta^lambda. If FALSE d~delta^lambda  
#'@param ... Further plot arguments passed: see 'plot.smacof' and 'plot' for detailed information.
#'
#'

#' 
#'@details
#' \itemize{
#' \item  Configuration plot (plot.type = "confplot"): Plots the MDS configuration.
#'  \item Residual plot (plot.type = "resplot"): Plots the dhats f(T(delta)) against the transformed fitted distances T(d(X)).
#'  \item (Linearized) Shepard diagram (plot.type = "Shepard"): Is shep.lin=TRUE a diagram with the transformed observed normalized dissimilarities (T(delta) on x)  against the transformed fitted distance (T(d(X) on y) as well as a loess curve and a regression line corresponding to type (linear without intercept for ratio, linear for interval and isotonic for ordinal). If shep.lin=FALSE it uses the untransformed delta. Note that the regression line corresponds to the optimal scaling results (dhat) only up to a linear transformation. 
#'  \item Transformation Plot (plot.type = "transplot"): Diagram with normalized observed dissimilarities (delta, light grey) and the normalized explicitly transformed dissimilarities (T(Delta), darker) against the untransformed fitted distances (d(X)) together with a nonlinear regression curve corresponding to the explicit transformation (fitted power transformation). This is most useful for ratio models with power transformations as the transformations can be read of directly. For other MDS models and stresses, it still gives a quick way to assess how the explicit transformations worked.  
#'  \item Stress decomposition plot (plot.type = "stressplot"): Plots the stress contribution in of each observation. Note that it rescales the stress-per-point (SPP) from the corresponding function to percentages (sum is 100). The higher the contribution, the worse the fit.
#'  \item Bubble plot (plot.type = "bubbleplot"): Combines the configuration plot with the point stress contribution. The larger the bubbles, the worse the fit.
#' \item histogram (‘plot.type = "histogram"’: gives a weighted histogram of the dissimilarities (weighted with tweightmat if exists else with weightmat). For optional arguments, see ‘wtd.hist’.
#' }
#'
#' @importFrom graphics plot text identify legend mtext par
#' @importFrom plotrix thigmophobe.labels
#' @importFrom stats loess lm predict coef
#' @importFrom grDevices chull
#'
#' @return no return value; just plots for class 'smacofP' (see details)
#' 
#' @export
#' 
#' @examples
#' dis<-as.matrix(smacof::kinshipdelta)
#' res<-powerStressMin(dis)
#' plot(res)
#' plot(res,"Shepard")
#' plot(res,"resplot")
#' plot(res,"transplot")
#' plot(res,"stressplot")
#' plot(res,"bubbleplot")
#' plot(res,"histogram")
plot.smacofP <- function (x, plot.type = "confplot", plot.dim = c(1, 2), bubscale = 1, col, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), hull.conf = list(hull = FALSE, col = 1, lwd = 1, ind = NULL), shepard.x=NULL, identify = FALSE, type = "p", cex=0.5, pch = 20, asp = 1, main, xlab, ylab, xlim, ylim, col.hist=NULL, legend = TRUE, legpos, loess=TRUE, shepard.lin=TRUE, ...)
{
     ## --- check type args:
    plot.type <- match.arg(plot.type, c("confplot", "Shepard", "transplot", "resplot","bubbleplot", "stressplot", "histogram"), several.ok = FALSE)

  ## --- check label lists
  if (is.null(label.conf$label)) label.conf$label <- TRUE
  if (is.null(label.conf$pos)) label.conf$pos <- 3
  if (is.null(label.conf$col)) label.conf$col <- 1
  if (is.null(label.conf$cex)) label.conf$cex <- 0.8
  if (identify) label.conf$label <- FALSE
  
  ## --- check hull list
  if (is.null(hull.conf$hull)) hull.conf$hull <- FALSE
  if (is.null(hull.conf$col)) hull.conf$col <- 1
  if (is.null(hull.conf$lwd)) hull.conf$lwd <- 1
  if (is.null(hull.conf$ind)) hull.conf$ind <- NULL
  if (hull.conf$hull && is.null(hull.conf$ind)) stop("Index vector for hulls needs to be specified!")

  ##------ Configuration plot 
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]

  #if (type == "n") label.conf$pos <- NULL
    if (plot.type == "confplot") {
        if (missing(col))  col <- 1
        if (missing(main)) main <- paste("Configuration Plot") else main <- main
        if (missing(xlab)) xlab <- paste("Dimension", x1, sep = " ") else xlab <- xlab
        if (missing(ylab)) ylab <- paste("Dimension", y1, sep = " ") else ylab <- ylab
        
        if (missing(xlim)) xlim <- range(x$conf[, x1])*1.1
        if (missing(ylim)) ylim <- range(x$conf[, y1])*1.1

        graphics::plot(x$conf[, x1], x$conf[, y1], main = main, type = type, 
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
            pch = pch, asp = asp, col = col, cex=cex, ...)

        if (label.conf$label){
            if(label.conf$pos==5){
                plotrix::thigmophobe.labels(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), cex = label.conf$cex, text.pos = NULL, col = label.conf$col)
                } else {
            graphics::text(x$conf[, x1], x$conf[, y1], labels = rownames(x$conf), cex = label.conf$cex, pos = label.conf$pos, col = label.conf$col)
                }
        }
        
        if (identify) {
            graphics::identify(x$conf[, x1], x$conf[, y1], labels = rownames(x$conf), cex = label.conf$cex, pos = label.conf$cex, col = label.conf$col)
        }

      if (hull.conf$hull) {
      ind <- hull.conf$ind
      n <- dim(x$conf)[1]
      M <- as.data.frame(x$conf)
      XX <- cbind(M, ind) 
      X.sort <- XX[order(ind), ]
      xx <- yy <- NULL
      k <- 0
      for (i in 1:n) { 
        v <- X.sort$ind[i+1]
        if (i==n) v="$"
        if (X.sort$ind[i] == v ) { 
          k<-k+1 
        } else { 
          von <- i-k 
          xx <- X.sort[von:i, 1] 
          yy <- X.sort[von:i, 2]
          hpts <- grDevices::chull(x = xx, y = yy)
          hpts <- c(hpts, hpts[1])
          lines(xx[hpts], yy[hpts], col = hull.conf$col, lwd = hull.conf$lwd) 
          k<-0 
        } 
       } 
     }
   }
  
  #---------------- Shepard diagram ------------------   
    if (plot.type == "Shepard") {

        if (missing(main)) {
            main <- paste("Shepard Diagram")
            if(shepard.lin) main <- paste("Linearized",main)
            } else main <- main
        if (missing(xlab)) {
            if (is.null(shepard.x)) xlab <- "Dissimilarities" else xlab <- "Proximities"
              if(shepard.lin) xlab <- paste("Transformed",xlab)
         } else xlab <- xlab

        if (missing(ylab)) ylab <- "Transformed Configuration Distances" else ylab <- ylab

        wm <- x$tweightmat
        if(is.null(wm)) wm <- x$weightmat 
        #additional optimal scaling to make the dhat usable
        #TODO: Eventually figure out how we can use the dhat$iord for r!=1 as they are the correct functions.  
        #type <- x$type
        #trans <- type
        #typo <- type
        #ties <- "primary"
        #if (trans=="ratio"){
        #trans <- "none"
        #}
        #else if (trans=="ordinal" & ties=="primary"){
        #trans <- "ordinalp"
        #typo <- "ordinal (primary)"
        #} else if(trans=="ordinal" & ties=="secondary"){
        #trans <- "ordinals"
        #typo <- "ordinal (secondary)"
        #} else if(trans=="ordinal" & ties=="tertiary"){
        #trans <- "ordinalt"
        #typo <- "ordinal (tertiary)"
        #}
        #disobj <- smacof::transPrep(as.dist(x$delta), trans = trans, spline.intKnots = 2, spline.degree = 2)
        #if(shepard.lin) disobj <- smacof::transPrep(as.dist(x$tdelta), trans = trans, spline.intKnots = 2, spline.degree = 2)
        #e <- as.dist(x$confdist)
        #wm <- x$tweightmat
        #if(is.null(wm)) wm <- x$weightmat
        #n <- x$nobj
        #dhat2 <- smacof::transform(e, disobj, w = as.dist(wm), normq = n )  ## dhat update
        #iord <- dhat2$iord.prim
        #dhatt <- dhat2$res
        #dhats <- structure(dhatt, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE)
       #FIXME: labels
       # dhats <- as.matrix(dhatd)
        notmiss <- as.vector(as.dist(x$weightmat) > 0)
        if (is.null(shepard.x)) {
           delts <- as.vector(x$delta) #with shepard.lin=FALSE we use the original delta
           if(shepard.lin) delts <- as.vector(x$tdelta) #dhat) #tdelta) #tdelta) #as.vector(x$tdelta) #with shepard.lin=FALSE we use the Shepard diagram on the level of the T(Delta) as we approx T(Delta) by the confdists 
          } else {
           delts <- as.vector(as.dist(shepard.x))
          }
        confd <- as.vector(x$confdist) #Confdist are already transformed
        wm <- as.vector(wm)
         #delts=xcoor in smacof 
         if (missing(xlim)) xlim <- range(delts[notmiss],na.rm=TRUE)
         if (missing(ylim)){
             ylim <- range(confd[notmiss])
             #ylim <- range(confd)
             ylim[1] <- 0
         }

        if(missing(col)) col <- c("grey70","grey40","black")
         
        #delta=observed delta Delta
        #tdelta=transformed delta normalized T(Delta) 
        #distances= dhats, optimally scaled transformed Delta and normalized f(T(Delta))
        graphics::plot(delts[notmiss], confd[notmiss], main = main, type = "p", pch=pch, cex = cex, xlab = xlab, ylab = ylab, col = col[1], xlim = xlim, ylim = ylim, ...)
        #notmiss.iord <- notmiss[x$iord]
        notmiss.iord <- notmiss[x$iord]
        delts1 <- delts[notmiss]
        confd1 <- confd[notmiss]
        wm1 <- wm[notmiss]
        dhats1 <- as.vector(x$dhat)[notmiss]
        expo <- 1
        disttrans.ind <- names(x$pars)%in%c("kappa","r") #TODO: make sure only the distance parameter is here, so kappa or r or whatever it is with the stops functions, enhance this with any new parameter names is it exists. ALso, if it is kappa or mu we take it at face vlaue and if it is r we need to double it as kappa=2*r.
        #TODO: Be careful not to name different parameter the same way as some of the distance transformation parameters
        disttrans <- x$pars[disttrans.ind]
        #points((delts[x$iord])[notmiss.iord], sqrt(2*x$nobj)*(as.vector(x$dhat[x$iord]))[notmiss.iord], type = "b", pch = pch, cex = cex,col=col[3])
        #TRIED: tried to change the sqrt(2*nobj) which works for r=0.5/kappa=1. I think part of the issue is that we also do enorm - so can we figure out a scaling factor from the scale of the confdist?
       #  SOLVED: used a linear model to get a scaling factor and an intercept. Looks good! The dhat are on a scale that is just transformed with enorm() and and we get the scaling factors with a linear function
        if(x$type=="ratio")
        {
            ## I checked with plot vs. confd~dhat: In case "ratio" this must go through 0 so we do not fit an intercept
        scallm <- stats::coef(stats::lm(confd1~-1+dhats1,weights=wm))
        scallm <- c(0,scallm)
        #dhatsscal <- as.vector(x$dhat[x$iord])
        #dhatsscal <- scallm[1]+scallm[2]*dhatsscal   
        }
        if(x$type=="interval")
        {
            ## In case of interval it needs not go through 0, so we fit intercept.
            scallm <- stats::coef(stats::lm(confd1~dhats1,weights=wm)) #with intercept is better for interval and ratio; and ordinal and kappa <=1. A bit less good at higher kappa, but overall better to have the intercept.
                                        #cat(scallm,"\n")
           #dhatsscal <- as.vector(x$dhat[x$iord])
           #dhatsscal <- scallm[1]+scallm[2]*dhatsscal
        }
        if(x$type=="ordinal")
        {
            #TODO: change here if we have objects with otehr parameters. First one must always be the configuration distance transformation (usually kappa, r, or mu) 
            expo <- switch(names(disttrans),
                       r=2*disttrans,
                       kappa=disttrans
                       )                   
            ## In case of ordinal I also use lm with intercept to get the scaling factor, but we also need to take the power transformation into account, so I create expo and do the lm thus.  
            ## If found that re-scaling of the dhats gets better if we do this power regression
            ## I'm not 100% sure why but I think it is because the isotonic regression is invariant to parametric transformations of the dhats which are out x argument in isoreg: so we need to manually include the power transformation somwhow. It may not be 100% correct because of enorm() but it looks better than ever. But in the the metric MDS the power transformation is taken into account by the tdelta.  
            scallm <- stats::coef(stats::lm(confd1~I(dhats1^expo),weights=wm))
            #scallm2 <- coef(lm(confd1^(1/expo)~dhats1,weights=wm)) #alternative where only the predictor gets transformed; looks a bit less but is a bit less accurate in trials. test more  
            #ir1 <- stats::isoreg(x=dhats1,y=confd1)
            #dhatscal+(ir1$yf[x$iord]-dhatsscal)    
        }
        if(x$type=="mspline")
        {
            #TODO: change here if we have objects with other parameters. First one must always be the configuration distance transformation (usually kappa, r, or mu) 
            expo <- switch(names(disttrans),
                       r=2*disttrans,
                       kappa=disttrans
                       )                    
            scallm <- stats::coef(stats::lm(confd1~I(dhats1^expo),weights=wm))
        }
        #scallm <- coef(lm(confd1~-1+dhats1,weights=wm))
        #scallm <- c(0,scallm)
        #cat(scallm,"\n")
        points((delts[x$iord])[notmiss.iord], scallm[1]+scallm[2]*(as.vector(x$dhat[x$iord])^expo)[notmiss.iord], type = "b", pch = pch, cex = cex,col=col[3])
        #Alternative: points((delts[x$iord])[notmiss.iord], (scallm2[1]+(scallm2[2]*(dhats1[x$iord])[notmiss.iord]))^expo, type = "b", pch = pch, cex = cex,col="red")
         #points((delts[x$iord])[notmiss.iord], (scallm1[1]+scallm1[2]*(dhats1[x$iord])[notmiss.iord])^expo, type = "b", pch = pch, cex = 2,col="green")#
        ##NOTE: I can't make smacofs transform work with normq=n in our fitting functions, so I scale up the dhat that are obtained from smacof::transform to the scale of the confdist that is returned.
        ## Since we we use normq=0.5 in fitting functions we thus need to scale the dhats up with sqrt(2*n)
        ## because in transform they do a=delta * sqrt(normq/sum(weights*delta^2)) and we want normq=n 
        ## so if we mutliply a*sqrt(2*n) it is as if we set normq=n.
        ## Still doesn't work because the r or kappa transformation isn't properly reflected and the x$dhats are only correct with k=1, r=0.5. It looks like there is some sort of scaling factor I'd have to apply but I don't know which one
        #delts1 <- delts[notmiss]
        #confd1 <- confd[notmiss]
        #wm1 <- wm[notmiss]
        if(loess) {
            ##no need to distinguish as in loess there is no constant 
            #if(x$type=="ratio") ptl <- predict(stats::loess(confd~-1+delts,weights=wm))
            #if(x$type=="interval") ptl <- predict(stats::loess(confd1~delts1,weights=wm))
            #if(x$type=="ordinal") ptl <- predict(stats::loess(confd1~delts1,weights=wm))
            ptl <- predict(stats::loess(confd1~delts1,weights=wm),family="symmetric")
            graphics::lines(delts[order(delts)],ptl[order(delts)],col=col[2],type="b",pch=pch,cex=cex)
        }
        #NOTE: This code would do the transformations manually based on f(confd~delts). That needs to coincide up to a scaling factor with the object$dhat, so I included this for checking that it works (mainly because the manual isoreg  and the isoreg in smacof do not give the same results and the former can't take weights), so I'd like to stick with the object$dhat as fitted in the MDS. For ratio and interval it would make no difference anyway.   
        #if(x$type=="ordinal")  {
        ## NOTE: we now do manual isotonic regression here as with our implementation the dhats from smacof are on a different scale. This is not 100% correct as we don't take the weightmat into account but for diagnostics its cool. 
        #   ir <- stats::isoreg(x=delts1,y=confd1) 
        #   #ptl <- ir$yf[ir$ord] 
        #  graphics::lines(ir,col=col[3],pch=pch,cex=cex,do.points=TRUE)
        #} else { 
        #if(x$type=="ratio") pt <- predict(stats::lm(confd1~-1+delts1,weights=wm))
        #if(x$type=="interval") pt <- predict(stats::lm(confd1~delts1,weights=wm))
        #graphics::lines(delts[order(delts)],pt[order(delts)],col=col[3],type="b",pch=pch,cex=cex,lwd=1)
        #}
    # Looks like we can't just use the dhat[iord] idea because normq is different in the calls, check that out too and alos because the scale of the confdist changes due to the power and the enorm. SOLVED: is there a relationhsip to figure out? Yes, linear no intercept for ratio, linear with intercept for interval, linear with pwoer function for ordinal       
    }
    if (plot.type == "transplot") {
        if(missing(col)) col <- c("grey40","grey70","grey30")#,"grey50")
        
        disttrans.ind <- names(x$pars)%in%c("kappa","r") #TODO: make sure only the distance parameter is here, so kappa or r or whatever it is with the stops functions, enhance this with any new parameter names is it exists. ALso, if it is kappa or mu we take oit ast face vlaue and if it is r we need to double it as kappa=2*r 
        disttrans <- x$pars[disttrans.ind]
        if(names(disttrans)%in%("r")) disttrans <- 2*disttrans 
        if(sum(disttrans.ind)==0)
             {
              warning("I can't identify the distance transformation parameter. I will use untransformed distances.")
              disttrans <- 1
             }
             deltao <- as.vector(x$delta/enorm(x$delta)) #normalize the delta
             deltat <- as.vector(x$tdelta/enorm(x$tdelta)) #are already normalized in smacofP but not in copsc. TODO: should we use the dhat here? so just explicit or both?
             dreal <- as.vector(x$confdist)^(1/disttrans) #change the confdist back to the eudclidean distance; we could also do dist(x$conf) but first is quicker 
             if (missing(main)) main <- paste("Transformation Plot")
             else main <- main
             if (missing(ylab)) ylab <- "Normalized Dissimilarities"
             else xlab <- xlab
             if (missing(xlab))  xlab <- "Untransformed Configuration Distances"
             else ylab <- ylab
            if (missing(ylim))  ylim <- c(min(min(deltat),min(deltao)),max(max(deltat),max(deltao)))
            #ylim[1] <- 0 
        if (missing(xlim)) xlim <- c(min(dreal),max(dreal))
            graphics::plot(dreal, deltao, main = main, type = "p", cex = 0.75, xlab = xlab, ylab = ylab, col = col[2], xlim = xlim, ylim = ylim, pch=20)
            #graphics::plot(deltat,dreal, main = main, type = "p", cex = 0.75, xlab = ylab, ylab = xlab, col = col[2], xlim = ylim, ylim = xlim,pch=20)
            graphics::points(dreal, deltat, type = "p", cex = 0.75, col = col[1],pch=20)
           #if(x$type=="ratio") pt <- predict(stats::lm(deltat~-1+I(dreal^disttrans))) #with intercept forcing thorugh 0
                                        #if(x$type=="interval"|| x$type=="ordinal") pt <- predict(stats::lm(deltat~I(dreal^disttrans))) #with intercept forcing through 0
            pt <- predict(stats::lm(deltat~I(dreal^disttrans))) #with intercept 
            #pt2 <- predict(stats::lm(deltat~I(dreal^kappa))) #with intercept not forcing thorugh 0 
            #po <- predict(stats::lm(deltao~-1+I(dreal^kappa))) #with intercept
            #lines(deltat[order(deltat)],pt[order(deltat)],col=col[1],type="b",pch=20,cex=0.5)
            #lines(deltao[order(deltao)],po[order(deltao)],col=col[2],type="b",pch=20,cex=0.5)
            #graphics::lines(dreal[order(dreal)],po[order(dreal)],col=col[4])
            graphics::lines(dreal[order(dreal)],pt[order(dreal)],col=col[3],type="b",pch=19,cex=0.1)
            #graphics::lines(dreal[order(dreal)],po[order(dreal)],col=col[4],type="b",pch=19,cex=0.25) 
            if(legend) {
                if(missing(legpos)) legpos <- "topleft" 
                graphics::legend(legpos,legend=c("Transformed","Untransformed"),col=col[1:2],pch=20)
            }
    }
 #--------------- Residual plot --------------------
    
if (plot.type == "resplot") {
        dhats <- as.vector(x$dhat)
        confd <- as.vector(x$confdist)
        
        if(missing(col)) col <- col <- c("grey40","grey70","black")
        if (missing(main)) main <- paste("Residual Plot") else main <- main
        if (missing(xlab)) xlab <- "Disparities (d-hats)" else xlab <- xlab
        if (missing(ylab)) ylab <- "Transformed Configuration Distances" else ylab <- ylab

        if (missing(xlim)) xlim <- range(c(0,dhats))
        if (missing(ylim)) ylim <- range(c(0,confd))
        graphics::plot(dhats, confd, main = main, 
            type = "p", col = col[1], xlab = xlab, ylab = ylab, 
            xlim = xlim, ylim = ylim, ...)
        abline(lm(confd~-1+dhats),col=col[3])
        if(loess) {
        ptl <- predict(stats::loess(confd~dhats))
        graphics::lines(dhats[order(dhats)],ptl[order(dhats)],col=col[2],type="b",pch=pch,cex=cex)
        }
#        abline(0,1,col="red")
}
#----------------------- Stress decomposition -----------------
        if (plot.type == "stressplot") {
        if(missing(col)) col <- "lightgray"
        if (missing(main)) main <- paste("Stress Decomposition Chart") else main <- main
        if (missing(xlab)) xlab <- "Objects" else xlab <- xlab
        if (missing(ylab)) ylab <- "Stress Proportion (%)" else ylab <- ylab

        spp.perc <- sort((x$spp/sum(x$spp) * 100), decreasing = TRUE)
        xaxlab <- names(spp.perc)

        if (missing(xlim)) xlim1 <- c(1, length(spp.perc)) else xlim1 <- xlim
        if (missing(ylim)) ylim1 <- range(spp.perc) else ylim1 <- ylim

        oldpar <- par(no.readonly=TRUE)
        on.exit(par(oldpar))
        par(mar = c(3, 4, 4, 2))
        plot(1:length(spp.perc), spp.perc, xaxt = "n", type = "p", 
            xlab = " ", ylab = ylab, main = main, xlim = xlim1, 
            ylim = ylim1, ...)
        graphics::mtext(xlab, side=1, padj=2)
        graphics::text(1:length(spp.perc), spp.perc, labels = xaxlab, pos = label.conf$pos, cex = 0.8)
        for (i in 1:length(spp.perc)) lines(c(i, i), c(spp.perc[i],0), col=col, lty = 2)
        }

    
  #------------------------------ bubble plot -------------------------
    if (plot.type == "bubbleplot") {
        if(missing(col)) col <- 1
        if (missing(main)) main <- paste("Bubble Plot") else main <- main
        if (missing(xlab)) xlab <- paste("Dimension", x1, sep = " ") else xlab <- xlab
        if (missing(ylab)) ylab <- paste("Dimension", y1, sep = " ") else ylab <- ylab

        if (missing(xlim)) xlim <- range(x$conf[, x1]) * 1.1
        if (missing(ylim)) ylim <- range(x$conf[, y1]) * 1.1

        spp.perc <- x$spp/sum(x$spp) * 100
        bubsize <- spp.perc/length(spp.perc) * (bubscale+3)

        graphics::plot(x$conf, cex = bubsize, main = main, xlab = xlab, 
            ylab = ylab, xlim = xlim, ylim = ylim, asp=asp, ...)
        xylabels <- x$conf
        ysigns <- sign(x$conf[, y1])
        xylabels[, 2] <- (abs(x$conf[, y1]) - (x$conf[, y1] * (bubsize/50))) * ysigns
        graphics::text(xylabels, rownames(x$conf), pos = label.conf$pos, cex = 0.7)
    }

  #------------------------------ histogram plot -------------------------
  if (plot.type == "histogram")
  {
    if (missing(main)) main <- paste("Weighted Histogram") else main <- main
    if (missing(xlab)) xlab <- paste("Dissimilarity") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Frequency") else ylab <- ylab

    wghts <- x$weightmat
    if(!is.null(x$tweightmat)) wghts <- x$tweightmat 
    weights::wtd.hist(x$delta, weight = wghts, main = main, xlab = xlab, ylab = ylab, col = col.hist, ...) 
  }
    
 }

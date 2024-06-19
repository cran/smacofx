
#'@importFrom graphics plot abline lines text identify legend points hist
#'@importFrom stats predict loess lm 
#'@export
plot.cmdscalex <- function(x, plot.type = "confplot", plot.dim = c(1, 2), col, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), hull.conf = list(hull = FALSE, col = 1, lwd = 1, ind = NULL), identify = FALSE, type = "p", cex=0.5, pch = 20, asp = 1, main, xlab, ylab, xlim, ylim, col.hist=NULL, legend=TRUE, legpos, loess=TRUE,...)
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
 
   ##------ Configuration plot  (like for smacofB)
    x1 <- plot.dim[1]
    y1 <- plot.dim[2]
    ## if (type == "n")
    ## label.conf$pos <- NULL
    
    if (plot.type == "confplot") {
        if(missing(col)) col <- 1
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
      if (plot.type %in% c("Shepard","resplot")) {
        delt <- as.vector(x$dhat) 
        confd <- as.vector(x$confdist) 
        if(missing(col)) col <- c("grey60","grey65") 
        if (missing(main)) 
            main <- ifelse(plot.type=="Shepard",paste("Linearized Shepard Diagram"),paste("Residual Plot"))
        else main <- main
        if (missing(xlab)) 
            xlab <- "Transformed Dissimilarities"
        else xlab <- xlab
        if (missing(ylab)) 
            ylab <- "Configuration Distances"
        else ylab <- ylab
        if (missing(xlim)) 
            xlim <- range(delt)
        if (missing(ylim)) 
            ylim <- range(confd)
        graphics::plot(delt, confd, main = main, 
            type = "p", pch = ifelse(plot.type=="Shepard",20,1), cex = ifelse(plot.type=="Shepard",0.75,1), xlab = xlab, ylab = ylab, 
            col = col[1], xlim = xlim, ylim = ylim, ...)
        if(plot.type=="Shepard") {
           pt <- predict(stats::loess(confd~delt))
           graphics::lines(delt[order(delt)],pt[order(delt)],col=col[2],type="b",pch=20,cex=0.5)
        }
     graphics::abline(stats::lm(confd~delt))
    }
    if (plot.type == "transplot") {
             if(missing(col)) col <- c("grey40","grey70","grey30","grey60")
             if(is.null(x$theta)) x$theta <- 1
             deltao <- as.vector(x$delta^(1/x$theta))
             deltat <- as.vector(x$delta)
             dreal <- as.vector(x$confdist)
             if (missing(main)) main <- paste("Transformation Plot")
             else main <- main
             if (missing(ylab)) ylab <- "Dissimilarities"
             else ylab <- ylab
             if (missing(xlab)) xlab <- "Configuration Distances"
             else xlab <- xlab
             if (missing(ylim)) ylim <- c(min(deltat,deltao),max(deltat,deltao))
             if (missing(xlim)) xlim <- range(as.vector(dreal))
            graphics::plot(dreal, deltao, main = main, type = "p", cex = 0.75, xlab = xlab, ylab = ylab, col = col[2], xlim = xlim, ylim = ylim, pch=20)
            graphics::points(dreal, deltat, type = "p", cex = 0.75, col = col[1],pch=20)
            pt <- predict(stats::lm(deltat~dreal))
            #po <- predict(stats::lm(deltao~dreal))
            graphics::lines(dreal[order(dreal)],pt[order(dreal)],col=col[3],type="b",pch=19,cex=0.1)
            #graphics::lines(dreal[order(dreal)],po[order(dreal)],col=col[4])
             if(legend)
                 {
                if(missing(legpos)) legpos <- "topleft" 
                graphics::legend(legpos,legend=c("Transformed","Untransformed"),col=col[1:2],pch=20)
                }
         }
    if (plot.type == "stressplot") {
        warning("Plot not supported for this class. Please use a stress-based MDS.")
    }
    if (plot.type == "bubbleplot") {
        warning("Plot not supported for this class. Please use a stress-based MDS.")
    }

  #------------------------------ histogram plot -------------------------
   if (plot.type == "histogram")
   {
    if (missing(main)) main <- paste("Weighted Histogram") else main <- main
    if (missing(xlab)) xlab <- paste("Dissimilarity") else xlab <- xlab
    if (missing(ylab)) ylab <- paste("Frequency") else ylab <- ylab
    graphics::hist(x$delta, main = main, xlab = xlab, ylab = ylab, col = col.hist, ...) 
  }
 }

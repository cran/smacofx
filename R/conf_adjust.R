





## #' A static 3d plot S3 generic 
## #'
## #' A static 3d plot 
## #'
## #' @title plot3dstatic: static 3D plots 
## #' @param x object
## #' @param plot.dim dimensions to plot
## #' @param main main title
## #' @param xlab label for x axis 
## #' @param ylab label for y axis
## #' @param zlab label for z axis
## #' @param col color
## #' @param ... other arguments
## #' 
## #' @export
## plot3dstatic <- function(x, plot.dim = c(1,2,3), main, xlab, ylab, zlab, col, ...) UseMethod("plot3dstatic")

## #' 3D plots: plot3dstatic method for class cmdscale
## #'
## #' 
## #' This methods produces a static 3D configuration plot.
## #' @param x object of class cmdscale
## #' @param plot.dim vector of length 3 with dimensions to be plotted
## #' @param main plot title
## #' @param xlab label of x axis
## #' @param ylab label of y axis
## #' @param zlab label of z axis
## #' @param col color of the text labels
## #' @param ... Further plot arguments passed: see 'scatterplot3d' in package 'scatterplot3d' for detailed information.
## #'
## #'@export
## #'@import scatterplot3d
## plot3dstatic.cmdscalex <- function (x, plot.dim = c(1, 2, 3), main, xlab, ylab, zlab, col,...) 
## {
##     ndim <- dim(x$points)[2]
##     if (ndim < 3) 
##         stop("No 3D plots can be drawn for ndim < 3 !")
##     if (length(plot.dim) != 3) 
##         stop("plot.dim must be of length 3!")
##     pd1 <- plot.dim[1]
##     pd2 <- plot.dim[2]
##     pd3 <- plot.dim[3]
##     if (pd3 > ndim) 
##         stop("Only", ndim, "dimensions were extracted!")
##     if (missing(xlab)) 
##         xlab <- paste("Dimension", pd1)
##     if (missing(ylab)) 
##         ylab <- paste("Dimension", pd2)
##     if (missing(zlab)) 
##         zlab <- paste("Dimension", pd3)
##     x1 <- x$points[, pd1]
##     y1 <- x$points[, pd2]
##     z1 <- x$points[, pd3]
##     if (missing(main)) 
##         main1 <- "Configuration Plot"
##     else main1 <- main
##     if (missing(col)) 
##         col1 <- "blue"
##     else col1 <- col
##     pr <- scatterplot3d::scatterplot3d(x1, y1, z1, type = "n", main = main1, xlab = xlab, ylab = ylab, zlab = zlab, ...)
##     text(pr$xyz.convert(x1, y1, z1), labels = rownames(x$points), col = col1)
## }

#'conf_adjust: a function to procrustes adjust two matrices
#'
#'@param conf1 reference configuration, a numeric matrix
#'@param conf2 another configuration, a numeric matrix 
#'@param verbose should adjustment be output; default to FALSE
#'@param eps numerical accuracy
#'@param itmax maximum number of iterations
#'@return a list with 'ref.conf' being the reference configuration, 'other.conf' the adjusted coniguration and 'comparison.conf' the comparison configuration
#'@export
conf_adjust<- function(conf1,conf2,verbose = FALSE,eps = 1e-12, itmax = 100)
 {
x0 <- conf1
n <- nrow(x0)
ndim <- 2
metric <- TRUE
xx <- conf2
kk <- diag(ndim)
cc <- matrix(0, n, ndim)
bb <- matrix(0, n, ndim)
yy <- xx
oloss <- Inf
itel <- 1
repeat {
y0 <- matrix(0, n, ndim)
y0 <- y0 + xx %*% kk
y0 <- ((n - 1) * y0)/(n * (n - 2))
zz <- matrix(0, n, ndim)
zz <- zz + xx %*% kk
xz <- crossprod(xx, zz)
kk <- procruster(xz)
nloss <- 0
 for (i in 1:n) {
   yy <- xx %*% kk    
   yy[i,] <- n * y0[i, ]/(n - 1)
   yy <- yy - outer(rep(1, n), y0[i, ]/(n - 1))
   nloss <- nloss + sum((y0 - yy)^2)
  }
        if (verbose) {
            cat("Iteration: ", formatC(itel, digits = 3, width = 3), 
                "Old Loss: ", formatC(oloss, digits = 10, width = 15, 
                  format = "f"), "New Loss: ", formatC(nloss, 
                  digits = 10, width = 15, format = "f"), "\n")
        }
        if (((oloss - nloss) < eps) || (itel == itmax)) {
            (break)()
        }
        itel <- itel + 1
        oloss <- nloss
}
x0 <- x0 %*% procruster(crossprod(x0,y0))
result <- list(ref.conf = x0, other.conf = yy, comparison.conf = y0)
result
}



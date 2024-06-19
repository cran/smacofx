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



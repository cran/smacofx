#' Corpse Paint
#'
#' A matrix of gray scale images of people in "corpse paint", a black-and-white make-up, plus a surprise.
#'
#' The images are gray scale 8 bit, i.e., 0-255 unique gray values scaled to be between 0 and 1. There are 32 total images with pixel size of 90 x 90 that have been vectorized to 32 columns labeled as "F1" through "F32". An image i can be reconstructed with matrix(corpsepaint[,i],ncol=90,nrow=90) 
#'
#' @format A 8100 x 32 matrix  
#' @name corpsepaint
#' @examples
#' oldpar<-par(no.readonly=TRUE)
#' par(mfrow=c(4,8))
#' for(i in 1:ncol(corpsepaint)){
#' p1<-matrix(corpsepaint[,i],ncol=90,nrow=90,byrow=FALSE) 
#' image(p1,col=gray.colors(256),main=colnames(corpsepaint)[i])
#' }
#' par(oldpar)
NULL

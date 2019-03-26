#' @useDynLib penda, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


#' Draw penda
#'
#' Take care, this function draw a panda, not penda, because it does not exists!
#'
#'@export
draw_penda = function(){
  tetey = c(4, 5, 5.8, 5.8, 5.4, 5, 4, 3, 1.5, 1.3, 1.5, 2.5, 4)
  tetex = c(1, 1.75, 2.9, 4, 4.8, 5.2, 5.8, 6, 5, 3.8, 2.3, 1.2, 1)
  orgy = c(4.5, 5.2, 6, 6.2, 5.8, 5, 4.5)
  orgx = c(1.35, 1, 1.2, 2.3, 2.9, 1.75, 1.35)
  ordy = c(5.8, 6.2, 6, 5.6, 4.8, 4.2, 5, 5.4, 5.8)
  ordx = c(4, 4.8, 5.3, 5.9, 6, 5.7, 5.2, 4.8, 4)
  
  ydx = c(3.9, 4.2, 4.8, 5, 4.8, 4.2, 3.8, 3.9)
  ydy = c(4, 4.4, 4.2, 4, 3.3, 3.2, 3.5, 4)
  
  ygx = c(2.5, 2.2, 2.5, 3.2, 3.5, 3.2, 2.5)
  ygy = c(3.2, 4, 4.3, 4.3, 4, 3.3, 3.2 )
  
  points = data.frame(x = c(2.8, 3.5, 4.3)
                      ,y = c(3.8, 2.8, 3.8))
  plot(0,0,col=0, 
    # main="PENDA: PErsoNalized Data Analysis"
    xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
    xlim=c(0,7), ylim=c(0,7)
  )
  # grid()
  polygon(tetex, tetey)
  polygon(ordx, ordy, lwd=8)
  polygon(orgx, orgy, lwd=8)
  polygon(ydx, ydy, lwd=8)
  polygon(ygx, ygy, lwd=8)
  points(points, cex=3, pch=16)    
}

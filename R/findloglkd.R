

# Function to find both log likelihood values of horizontal and vertical gradient
# given observed blurred image

findloglkd <- function(y, kern = c("norm","circnorm","cauchy","disc"),
                    h = seq(1,5.5,by = 0.05), sigma = seq(0.01,0.2,by = 0.01), ar = FALSE) {

  param <- expand.grid(h = h,sigma = sigma)

  dx.h <- rip.conv(y, rip.grad$x,"valid")
  dx.v <- rip.conv(y, rip.grad$y,"valid")
  Xh <- Mod(rip.ndft(dx.h))
  Xv <- Mod(rip.ndft(dx.v))

  if(ar) {
    G <- g.autoreg(y)
    Hh <- h.theoretical(dim(Xh))$h
    Hv <- h.theoretical(dim(Xv))$v
    likval <- apply(param,1,FUN = function(param){
      c(lkd_gen(Xh,blurkernel(kern,rad = param[1]),G = G$h,H = Hh,sigma = param[2]),
        lkd_gen(Xv,blurkernel(kern,rad = param[1]),G = G$v,H = Hv,sigma = param[2]))  } )
  } else {
    likval <- apply(param,1,FUN = function(param){
      c(lkd_gen(Xh,blurkernel(kern,rad = param[1]),sigma = param[2]),
        lkd_gen(Xv,blurkernel(kern,rad = param[1]),sigma = param[2])) } )
  }

  likdata <- data.frame(param, likh = likval[1,],likv = likval[2,])

  return(likdata)
}

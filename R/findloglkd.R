

# Function to calculate log likelihood for both horizontal and vertical gradient given observed blurred image
# over a sqaure rectangular of value for prior sd sigma and prior parameter
# Parameters:
# kern: "norm" for rectangular gaussian kernel,"circnorm" for circular gaussian kernel,
#       "cauchy" for circular cauchy kernel and "disc" for disc kernel
# rad: sequence of radius for circular support for circular kernels
# h: sequence of scale for blur kernel
# ar: If TRUE uses AR Prior

findloglkd <- function(y, kern = c("norm","circnorm","cauchy","disc"),rad = seq(1,5.5,by = 0.05),
                         h = rad, sigma = seq(0.01,0.2,by = 0.01), kap = 1, ar = FALSE, eta = 0.001) {

  dx.h <- rip.conv(y, rip.grad$x,"valid")
  dx.v <- rip.conv(y, rip.grad$y,"valid")
  Xh <- Mod(rip.ndft(dx.h))
  Xv <- Mod(rip.ndft(dx.v))

  if(ar){

    G <- g.autoreg(y)
    Hh <- h.theoretical(dim(Xh))$h
    Hv <- h.theoretical(dim(Xv))$v

    if(!kern == "norm"){  #If Circular support type Kernel

      param <- expand.grid(rad = rad,sigma = sigma)

      likval <- apply(param,1,FUN = function(param){
        c(lkd_gen(Xh,blurkernel(kern,rad = param[1],kap = kap),G = G$h,H = Hh,sigma = param[2],eta = eta),
          lkd_gen(Xv,blurkernel(kern,rad = param[1],kap = kap),G = G$v,H = Hv,sigma = param[2],eta = eta))})

    }else{   #If rectangular support type Kernel

      param <- expand.grid(h = h,sigma = sigma)

      likval <- apply(param,1,FUN = function(param){
        c(lkd_gen(Xh,blurkernel(kern,h = param[1]),G = G$h,H = Hh,sigma = param[2],eta = eta),
          lkd_gen(Xv,blurkernel(kern,h = param[1]),G = G$v,H = Hv,sigma = param[2],eta = eta))})
    }

   }else {

     if(!kern == "norm"){  #If Circular support type Kernel

       param <- expand.grid(rad = rad,sigma = sigma)

       likval <- apply(param,1,FUN = function(param){
         c(lkd_gen(Xh,blurkernel(kern,rad = param[1],kap = kap),sigma = param[2],eta = eta),
           lkd_gen(Xv,blurkernel(kern,rad = param[1],kap = kap),sigma = param[2],eta = eta))})

     }else{   #If rectangular support type Kernel

       param <- expand.grid(h = h,sigma = sigma)

       likval <- apply(param,1,FUN = function(param){
         c(lkd_gen(Xh,blurkernel(kern,h = param[1]),sigma = param[2],eta = eta),
           lkd_gen(Xv,blurkernel(kern,h = param[1]),sigma = param[2],eta = eta))})
     }

  }

  likdata <- data.frame(param, likh = likval[1,],likv = likval[2,])

  return(likdata)
}

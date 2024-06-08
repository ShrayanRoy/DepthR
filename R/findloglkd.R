

# Function to calculate log likelihood for both horizontal and vertical gradient given observed blurred image
# over a sqaure rectangular of value for prior sd sigma and prior parameter
# Parameters:
# kern: "norm" for Rectangular Gaussian kernel,"circnorm" for Circular Gaussian kernel,
#       "cauchy" for Circular Cauchy kernel and "disc" for Disc kernel
# rad: sequence of radius for circular support for circular kernels
# h: sequence of scale for blur kernel
# ar: If TRUE uses AR Prior

findloglkd <- function(y, kern = c("norm", "circnorm", "cauchy", "disc", "tcauchy"),rad = seq(1,5.5,by = 0.05),h = NULL,
                       sigma = seq(0.01,0.2,by = 0.01), kap = 1, ar = FALSE, eta = 0.010,thres = 0) {

  dx.h <- rip.conv(y, rip.grad$x,"valid")
  dx.v <- rip.conv(y, rip.grad$y,"valid")
  Xh <- Mod(rip.ndft(dx.h))
  Xv <- Mod(rip.ndft(dx.v))

  if(ar){

    G <- g.autoreg(y)
    Hh <- h.theoretical(dim(Xh))$h
    Hv <- h.theoretical(dim(Xv))$v

    if(kern != "norm"){  #If Circular support type Kernel

      param <- expand.grid(rad = rad,sigma = sigma)

      likval <- apply(param,1,FUN = function(param){
        c(lkd_gen(Xh,blurkernel(kern,rad = param[1],kap = kap),G = G$h,H = Hh,sigma = param[2],eta = eta,thres = thres),
          lkd_gen(Xv,blurkernel(kern,rad = param[1],kap = kap),G = G$v,H = Hv,sigma = param[2],eta = eta,thres = thres))})

    }else{   #If rectangular support type Kernel

      param <- expand.grid(h = h,sigma = sigma)

      likval <- apply(param,1,FUN = function(param){
        c(lkd_gen(Xh,blurkernel(kern,h = param[1]),G = G$h,H = Hh,sigma = param[2],eta = eta,thres = thres),
          lkd_gen(Xv,blurkernel(kern,h = param[1]),G = G$v,H = Hv,sigma = param[2],eta = eta,thres = thres))})
    }

   }else {

     if(kern != "norm"){  #If Circular support type Kernel

       param <- expand.grid(rad = rad,sigma = sigma)

       likval <- apply(param,1,FUN = function(param){
         c(lkd_gen(Xh,blurkernel(kern,rad = param[1],kap = kap),sigma = param[2],eta = eta,thres = thres),
           lkd_gen(Xv,blurkernel(kern,rad = param[1],kap = kap),sigma = param[2],eta = eta,thres = thres))})

     }else{   #If rectangular support type Kernel

       param <- expand.grid(h = h,sigma = sigma)

       likval <- apply(param,1,FUN = function(param){
         c(lkd_gen(Xh,blurkernel(kern,h = param[1]),sigma = param[2],eta = eta,thres = thres),
           lkd_gen(Xv,blurkernel(kern,h = param[1]),sigma = param[2],eta = eta,thres = thres))})
     }

  }

  likdata <- data.frame(param, likh = likval[1,],likv = likval[2,])

  return(likdata)
}


# Function to estimate blur kernel parameters based on continuous grid search
# yimg: Input observed image
# kern: Used blur kernel model
# interval: Search Interval for radius

BlurParamEst <- function(yimg,kern = c("norm", "circnorm", "cauchy", "disc", "tcauchy"),kap = 1,
                         rho = c(0.3,0.6), interval = c(1,8), sigma = 0.10, eta = 0.005,thres = 1e-05){

  dx.h <- rip.conv(yimg, rip.grad$x,"valid")
  dx.v <- rip.conv(yimg, rip.grad$y,"valid")
  Xh <- Mod(rip.ndft(dx.h))
  Xv <- Mod(rip.ndft(dx.v))
  G <- g.autoreg(yimg,rho = rho)
  Hh <- h.theoretical(dim(Xh))$h
  Hv <- h.theoretical(dim(Xv))$v

  lkdfunc <- function(param){
    kk <- blurkernel("cauchy",rad = param,kap =  1)
    likval <- lkd_gen(Xh,kk = kk,G = G$h,H = Hh,sigma = sigma,eta = eta,thres = thres) +
      lkd_gen(Xv,kk = kk,G = G$v,H = Hv,sigma = sigma,eta = eta,thres = thres)

    return(likval)
  }

  gridStepSearch(lkdfunc,lower_bound = interval[1],upper_bound = interval[2],max_iter = 10)$x

}

# Optimal Choice of parameters obtained from experiments
# Less noisy image: sigma = 0.10, eta = 0.005
# More noisy image: sigma = 0.10, eta = 0.010






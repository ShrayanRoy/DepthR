


# Utilities for finding blur kernel based on idea of decorrelation penalty
# The idea is adhoc, more theoretical development is needed


# Function for calculation of deconvolution step of the algorithm

deconv1_step <- function(x,k,G = 1,H = 1,sigma = 0.10,
                        eta = 0.010, thres = 1e-05){
  if(nrow(x)%%2 == 0){
    x <- x[-1,]
    if(!identical(G,1))  G <- G[-1,]
    if(!identical(H,1))  H <- H[-1,]
  }

  if(ncol(x)%%2 == 0){
    x <- x[,-1]
    if(!identical(G,1))  G <- G[,-1]
    if(!identical(H,1))  H <- H[,-1]
  }

  K <- rip.dft(rip.pad(k,pad = dim(x)))
  Denom <- ((eta^2)*H + (sigma^2)*G*(Mod(K)^2))
  Denom[Denom < thres] <- thres
  Fac <- ((sigma^2)*K)/Denom

  z <- rip.dft(rip.dft(x)*Fac,inverse = TRUE)

  return(rip.shift(z))
}

# Function to estimate blur kernel based on observed blurred image

find_decorr <- function(yblur,kern = c("norm","circnorm","cauchy","tcauchy","disc"),radseq = seq(1,7,by = 0.10),
                   kap = 1,lag.max = 20,rho = c(0.3,0.3),sigma = 0.10,eta = 0.010,thresh = 1e-05){

  xh <- rip.conv(yblur,k = rip.grad$x,"valid")
  xv <- rip.conv(yblur,k = rip.grad$y,"valid")

  G <- g.autoreg(yblur,rho = rho)
  Hh <- h.theoretical(dim(xh))$h
  Hv <- h.theoretical(dim(xv))$v

  pena <- NULL

  for(i in 1:length(radseq)){

    kk <- blurkernel(kern = kern,rad = radseq[i],kap = kap)

    zh <- deconv1_step(x = xh,k = kk,G = G$h,H = Hh,sigma = sigma,eta = eta,thres = 1e-05)
    zv <- deconv1_step(x = xv,k = kk,G = G$v,H = Hv,sigma = sigma,eta = eta,thres = 1e-05)

    rh.along <- acf(as.vector(t(zh)),lag.max = lag.max,plot = F)
    rh.across <- acf(as.vector(zh),lag.max = lag.max,plot = F)
    rv.along <- acf(as.vector(zv),lag.max = lag.max,plot = F)
    rv.across <- acf(as.vector(t(zv)),lag.max = lag.max,plot = F)

    rh.along <- rh.along$acf[-1]
    rh.across <- rh.across$acf[-1]
    rv.along <- rv.along$acf[-1]
    rv.across <- rv.across$acf[-1]

    pena <- c(pena, sum(rh.along^2) + sum(rv.along^2) + sum(rh.across^2) + sum(rv.across^2))
  }

  return(data.frame(radseq = radseq, pena = pena))
}




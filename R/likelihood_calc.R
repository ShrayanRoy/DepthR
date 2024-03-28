

# Function to calculate log likelihood value given fourier transform and kernel
# X : Discrete Fourier Transform of image gradient
# kk : Input kernel
# sigma : prior parameter
# eta : noise parameter
# thres: for numerical stability

lkd_gen <- function(X,kk,G = 1,H = 1,sigma = 0.05,eta = 0.001,thres = 0){

  if(nrow(X)%%2 == 0){
    X <- X[-1,]
    if(!identical(G,1))  G <- G[-1,]
    if(!identical(H,1))  H <- H[-1,]
  }

  if(ncol(X)%%2 == 0){
    X <- X[,-1]
    if(!identical(G,1))  G <- G[,-1]
    if(!identical(H,1))  H <- H[,-1]
  }

  kkpad <- rip.pad(kk,pad = dim(X))
  K <- Mod(rip.dft(kkpad))

  varX <- as.vector(((sigma^2)*(K^2)*G) + ((eta^2)*H))

  # adjustment for dividing by very small number
  varX[varX < thres] <- thres
  l <- dexp(as.vector(X^2),rate = 1/varX,log = TRUE)

  return(sum(l))
}

# Functions to calculate h_{\omega}
# Imported from https://github.com/deepayan/rip/blob/main/rip.recover/R/freqvar.R

h.1d <- function(N, omega = seq_len(N) - 1)
{
  Mod(complex(modulus = 1, argument= -2 * pi * omega / N) - 1)
}

h.theoretical <- function(d)
{
  list(h = as.rip(outer(rep(1, d[1]), h.1d(d[2]))),
       v = as.rip(outer(h.1d(d[1]), rep(1, d[2]))))
}


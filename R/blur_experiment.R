

# Function to select random patch from given image y
randPatch <- function(y,psize = 51){

  L <- dim(y)
  f1 <- sample(1:(L[1] - psize + 1),size = 1)
  f2 <- sample(1:(L[2] - psize + 1),size = 1)
  ypatch <- as.rip(y[f1:(f1 + psize - 1),f2:(f2 + psize - 1)])

  return(ypatch)
}

# Function to plot radhat vs. sigma
plotExp <- function(ypatch,kern = c("norm", "circnorm", "cauchy", "disc", "tcauchy"),rad = 3,
                    kap = 1, radseq = seq(0.5,6,by = 0.05),sigmaseq = seq(0.01,0.4,by = 0.02)){

  # Currently kern = "norm" not supported

  yblur <- rip.conv(k = blurkernel(kern = kern,rad = rad,kap = kap),ypatch/255,"valid")
  res <- findloglkd(yblur,kern = kern,ar = TRUE,rad = radseq,sigma = sigmaseq,kap = kap)

  radhat <- function(d,s){
    with(subset(d,sigma == s),rad[which.max(likh + likv)])
  }

  rad_est <- sapply(sigmaseq,FUN = radhat,d = res)

  plot(rad_est ~ sigmaseq,type = "o",xlab = "Sigma",ylab = "Estimated rad",pch = 19,
       col = "cornflowerblue",lwd = 2,main = paste("True rad = ",rad," & kernel = ",kern))
  abline(h = rad,lty = 2,col = "black")

}


#doblurExpr <- function(y,psize = 51,kern = c("norm","circnorm","cauchy","disc"),rad = 3,kap = 1,
#                       radseq = seq(1,6,by = 0.05),sigmaseq = seq(0.01,0.4,by = 0.025)){

  # Currently kern = "norm" not supported
#  L <- dim(y)
#  f1 <- sample(1:(L[1] - psize + 1),size = 1)
#  f2 <- sample(1:(L[2] - psize + 1),size = 1)
#  ypatch <- y[f1:(f1 + psize - 1),f2:(f2 + psize - 1)]

#  yblur <- rip.conv(k = blurkernel(kern = kern,rad = rad,kap = kap),ypatch,"valid")
#  res <- findloglkd(yblur,kern = kern,ar = TRUE,rad = radseq,sigma = sigmaseq,kap = kap)

#  radhat <- function(d,s){
#    with(subset(d,sigma == s),rad[which.max(likh + likv)])
#  }

#  rad_est <- sapply(sigmaseq,FUN = radhat,d = res)

#  plot(rad_est ~ sigmaseq,type = "o",xlab = "Sigma",ylab = "Estimated rad",pch = 19,
#      col = "cornflowerblue",lwd = 2,main = paste("True rad = ",rad," & kernel = ",kern))
#  abline(h = rad,lty = 2,col = "black")

#}

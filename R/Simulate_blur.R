
# Function to simulate defocus blur given input image and blurmap
# yimg: Input Image
# blurMap: blurMap i.e. matrix of rad values
# kern: Choice of kernel
# kap: same as in blurkernel

Simulate_blur <- function(yimg,blurMap,kern = c("norm", "circnorm", "cauchy", "disc", "tcauchy"),kap = 1){

  kern <- match.arg(kern)
  blurmap <- round(blurMap,2) # to reduce time

  radMap <- unique(as.vector(blurMap))
  yblurred <- yimg

  if(length(dim(yimg)) == 3){
    ## For RGB Channel Image
    for(i in 1:length(radMap)){
      kk <- blurkernel(kern = kern, rad = radMap[i],kap = kap)
      temp1 <- rip.conv(yimg[,,1],kk,"same")
      temp2 <- rip.conv(yimg[,,2],kk,"same")
      temp3 <- rip.conv(yimg[,,3],kk,"same")

      mask <- (blurMap == radMap[i])
      yblurred[,,1][mask] <- temp1[mask]
      yblurred[,,2][mask] <- temp2[mask]
      yblurred[,,3][mask] <- temp3[mask]
    }

  }else{
    ## For BW Image
    for(i in 1:length(radMap)){
      kk <- blurkernel(kern = kern,rad = radMap[i],kap = kap)
      temp <- rip.conv(yimg,kk,"same")

      mask <- (blurMap == radMap[i])
      yblurred[mask] <- temp[mask]
    }

  }

  return(as.rip(yblurred))

}






# Function to find Cropped Image Segment based on Segment Anything Model (SAM) output
# y: original grayscale image in matrix/rip format
# box: same as `bbox` in SAM output
# mask: same as `segmentation` in SAM output
# seg: logical flag. If TRUE then returns cropped image segment based on both box and mask and
#       If FALSE only cropped image based on box

# See https://github.com/facebookresearch/segment-anything/blob/main/segment_anything/automatic_mask_generator.py
# for details related to `box` and `mask` argument

findSeg <- function(y,box,mask = NULL,seg = T){

  # To adjust for indexes starting from 0 in python
  ycrop <- y[(box[[2]] + 1):(box[[2]] + box[[4]] + 1),
             (box[[1]] + 1):(box[[1]] + box[[3]] + 1)]

  if(seg){
    maskcrop <- mask[box[[2]]:(box[[2]] + box[[4]]),
                     box[[1]]:(box[[1]] + box[[3]])]
    ycrop[!maskcrop] <- 0
  }

  return(as.rip(ycrop))
}

# Function to Spatially Deblur an Image based on blur map
# img: Original RGB channel image
# blurMap: A matrix of blurkernel of parameters similar in size to `img` dimension
# lamb: Similar to rip.deconv `lambda` argument

# For more details see https://github.com/deepayan/rip/blob/main/rip.recover/R/nonblind.R

spatDeblur <- function(img, blurMap, kern = c("norm","circnorm","cauchy","disc"),
                       kap = 1, lamb = 0.01) {

  radMap <- unique(as.vector(blurMap))
  ydeblur <- img

  for(i in 1:length(radMap)) {
    kk <- blurkernel(kern = kern, rad = radMap[i],kap = kap)
    temp1 <- rip.deconv(as.rip(img[,,1]),k = kk, method = "direct", lambda = lamb,
                        rho = list(along = 0, across = 0), patch = 150, verbose = TRUE)
    temp1[temp1 > 1] <- 1

    temp2 <- rip.deconv(as.rip(img[,,2]),k = kk, method = "direct", lambda = lamb,
                        rho = list(along = 0, across = 0), patch = 150, verbose = TRUE)
    temp2[temp2 > 1] <- 1

    temp3 <- rip.deconv(as.rip(img[,,3]),k = kk, method = "direct", lambda = lamb,
                        rho = list(along = 0, across = 0), patch = 150, verbose = TRUE)
    temp3[temp3 > 1] <- 1

    mask <- (blurMap == radMap[i])
    ydeblur[,,1][mask] <- temp1[mask]
    ydeblur[,,2][mask] <- temp2[mask]
    ydeblur[,,3][mask] <- temp3[mask]
  }
  return(ydeblur)
}






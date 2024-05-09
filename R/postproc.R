

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
    maskcrop <- mask[box[[2]]:(box[[2]] + box[[4]] + 1),
                     box[[1]]:(box[[1]] + box[[3]] + 1)]
    ycrop[!maskcrop] <- 0
  }

  return(as.rip(ycrop))
}

# Function to Spatially Deblur RGB Image or Gray Scale Image based on blur map
# img: Original RGB channel Image or Gray Scale Image as rip object
# blurMap: A matrix of blur kernel of parameters similar in size to `img` dimension
# lamb: Similar to rip.deconv `lambda` argument
# alpha: The hyper-Laplacian parameter. alpha = 2 corresponds to Gaussian, and alpha = 1 to double
#  exponential or Laplacian. alpha = 0.8 is commonly used to model natural images, and often
#  referred to as a "sparse" prior because it puts relatively more weight to 0 gradients.

# For more details see https://github.com/deepayan/rip/blob/main/rip.recover/R/nonblind.R

spatDeblur <- function(img, blurMap, kern = c("norm", "circnorm", "cauchy", "disc", "tcauchy"),
                        kap = 1, lamb = 0.01,alpha = 2) {

  radMap <- unique(as.vector(blurMap))
  ydeblur <- as.array(img)

  for(i in 1:length(radMap)) {
    kk <- blurkernel(kern = kern, rad = radMap[i],kap = kap)
    temp <- rip.deconv(as.rip(img),k = kk, method = "direct", lambda = lamb, patch = 150,
                       rho = list(along = 0, across = 0), verbose = TRUE,alpha = alpha)

    mask <- (blurMap == radMap[i])
    ydeblur[mask] <- as.array(temp)[mask]
  }
  return(as.rip(ydeblur))
}

# Function to Spatially Deblur RGB Image or Gray Scale Image based on blur map
# Using Richardson-Lucy Algorithm
# img: Original RGB channel Image or Gray Scale Image as rip object
# blurMap: A matrix of blur kernel of parameters similar in size to `img` dimension
# lamb: Similar to rip.deconvlucy `niter` argument

# For more details see https://github.com/deepayan/rip/blob/main/rip.recover/R/nonblind.R

spatDeblurLucy <- function(img, blurMap, kern = c("norm", "circnorm", "cauchy", "disc", "tcauchy"),
                        kap = 1, niter = 25) {

  radMap <- unique(as.vector(blurMap))
  ydeblur <- as.array(img)

  for(i in 1:length(radMap)) {
    kk <- blurkernel(kern = kern, rad = radMap[i],kap = kap)
    temp <- rip.deconvlucy(as.rip(img),k = kk,niter = niter)

    mask <- (blurMap == radMap[i])
    ydeblur[mask] <- as.array(temp)[mask]
  }
  return(as.rip(ydeblur))
}






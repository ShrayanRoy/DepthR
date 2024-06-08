

# Utilities for finding blur kernel based on idea of decorrelation loss
# In the case of spatially varying blur

# Function for calculation of deconvolution step of the algorithm
# x: gradient of observed blurred image
# k: blur kernel
# sigma : prior parameter
# eta : noise parameter
# thres: for numerical stability

deconv2_step <- function(x,k,G = 1,H = 1,sigma = 0.10,
                        eta = 0.01,thres = 1e-05){
  K <- rip.dft(rip.pad(k,pad = dim(x)))

  Denom <- ((eta^2)*H + (sigma^2)*G*(Mod(K)^2))
  Denom[Denom < thres] <- thres
  Fac <- ((sigma^2) * K * G) / Denom

  z <- rip.dft((rip.dft(x) * Fac), inverse = TRUE)

  return(rip.shift(z))
}

# Function to estimate blur radius map based on observed blurred image
# based on Decorrelation loss
# yimg: Input observed blurred image
# boxdata: List of `bbox` for each segment identified by SAM
# maskdata: List of `segmentation` for each segment identified by SAM
# kern: Choice of kernel model
# radseq: sequence of radius for circular support for circular kernels
# lag.max: Maximum lag to consider in calculation of decorrelation loss
# rho: same as rho in AR prior

FindBlurMap_decorr <- function(yimg, boxdata, maskdata,
                        kern = c("norm","circnorm","cauchy","tcauchy","disc"),
                        radseq = seq(1, 8, by = 0.1), kap = 1,
                        lag.max = 20, sigma = 0.10, eta = 0.01, rho = c(0.3, 0.6))
{

  kern <- match.arg(kern)
  ## "Same" to avoid confusion using SAM output
  dx.h <- rip.conv(yimg, rip.grad$x,"same")
  dx.v <- rip.conv(yimg, rip.grad$y,"same")

  rmap <- NULL
  rloss <- list()
  ## par(mfrow = c(4,5))
  for(i in 1:length(boxdata)){

    cat("\rsegment ", i , " of ", length(boxdata))
    ## Bounding box for Segment
    dx.h_part <- findSeg(dx.h, box = boxdata[[i]], seg = FALSE)
    dx.v_part <- findSeg(dx.v, box = boxdata[[i]], seg = FALSE)

    ## For AR Prior
    Gh <- g.autoreg(cbind(dx.h_part, 0), rho = rho)$h
    Gv <- g.autoreg(rbind(dx.v_part, 0), rho = rho)$v
    Hh <- h.theoretical(dim(dx.h_part))$h
    Hv <- h.theoretical(dim(dx.v_part))$v

    pena <- numeric(length(radseq))
    ## mask
    mask <- findSeg(maskdata[[i]], box = boxdata[[i]], seg = FALSE)

    for (j in 1:length(radseq)){

      kk = blurkernel(kern, rad = radseq[j], kap = kap)
      z.h <- deconv2_step(dx.h_part, k = kk, G = Gh, H = Hh, sigma = sigma, eta = eta)
      z.v <- deconv2_step(dx.v_part, k = kk, G = Gv, H = Hv, sigma = sigma, eta = eta)

      ## To select elements in direction of gradient
      dx.along <- t(z.h)[t(mask)]
      dy.along <- z.v[mask]
      dx.across <- t(z.v)[t(mask)]
      dy.across <- z.h[mask]

      rx.along <- acf(dx.along, lag.max = lag.max, plot = FALSE)$acf[-1]
      ry.along <- acf(dy.along, lag.max = lag.max, plot = FALSE)$acf[-1]
      rx.across <- acf(dx.across, lag.max = lag.max, plot = FALSE)$acf[-1]
      ry.across <- acf(dy.across, lag.max = lag.max, plot = FALSE)$acf[-1]

      pena[j] <- sum(rx.along^2) + sum(ry.along^2) + sum(rx.across^2) + sum(ry.across^2)
    }

    rloss[[i]] <- data.frame(segment = i, radius = radseq, loss = pena)
    rmap[i] <- radseq[which.min(pena)]
  }
  cat("\n")

  blurMap = matrix(min(rmap),nrow = dim(yimg)[1],ncol = dim(yimg)[2])
  for(j in length(boxdata):1){
    blurMap[maskdata[[j]]] = rmap[j]
  }

  return(list(rmap = rmap,
              rloss = do.call(rbind, rloss),
              blurmap = blurMap))
}

# Function to find  coloured red map based on estimated radius map
# blurmap: Estimated blur radius map
# limits: controls the range of scale in plotted blurmap

Draw_ColourMap <- function(blurmap,limits = c(1,8)){

  graph_df <- data.frame(expand.grid(x = 1:(nrow(blurmap)),y = 1:(ncol(blurmap))),
                         radmap = c(blurmap))

  return(ggplot(graph_df, aes(x = y, y = nrow(blurmap) - x + 1)) +
           geom_tile(aes(fill = radmap)) + labs(x = " ",y = " ") +
           scale_fill_distiller(palette = "YlGnBu",limits = limits) + theme_classic())
}

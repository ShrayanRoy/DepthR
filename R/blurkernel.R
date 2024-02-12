blurkernel <- function(kern = c("norm","circnorm","cauchy","disc"),
                       rad = 5,h = rad,gridsize = 10){

  kdim <- max(floor(h - 0.5) + 1.5,0.5)
  if(!kern == "norm") kdim <- max(floor(rad - 0.5) + 1.5,0.5)

  x <- seq(-kdim,kdim,length.out = 2*kdim*gridsize)
  f <- cut(x,seq(-kdim,kdim,by = 1),include.lowest = TRUE)
  g <- expand.grid(f1 = f, f2 = f, KEEP.OUT.ATTRS = FALSE)

  if(kern == "norm"){
    g$kval <- as.vector(outer(dnorm(x,sd = h/3),dnorm(x,sd = h/3)))
  }else if(kern == "circnorm"){
    g$kval <- as.vector(outer(dnorm(x,sd = h/3),dnorm(x,sd = h/3)) * (outer(x^2,x^2,"+") < rad^2))
  }else if(kern == "cauchy"){
    g$kval <- as.vector((1/(outer(x^2,x^2,"+") + (h/3)^2)^1.5) * (outer(x^2,x^2,"+") < rad^2))
  }else if(kern == "disc"){
    g$kval <- as.vector(outer(x^2,x^2,"+") <= rad^2)
  }else{
    stop("Invalid Kernel Type Specified")
  }

  kk <- xtabs(kval ~ f1 + f2, g) / xtabs( ~ f1 + f2, g)
  kk[] <- kk/sum(kk)

  as.rip(kk)
}

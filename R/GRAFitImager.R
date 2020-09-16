GRAFitImager <- function(image = NULL, InitialModel = NULL, FinalModel = NULL, ModelNoise = NULL, Data = NULL) {

  cmap = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
  medimg = median(abs(image[image > 0 & Data$region]))
  maximg = max(abs(image[Data$region]))
  stretchscale = 1/medimg
  par(mfrow=c(1,4))
  magimage(image, stretchscale = stretchscale, lo = -maximg, 
           hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
           xlab = "x/pix", ylab = "y/pix")
  mtext("Data", side=3) 
  
  magimage(InitialModel, stretchscale = stretchscale, lo = -maximg, 
           hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
           xlab = "x/pix", ylab = "y/pix")
  mtext("Initial Model", side=3)  
  
  magimage(FinalModel, stretchscale = stretchscale, lo = -maximg, 
           hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
           xlab = "x/pix", ylab = "y/pix")
  mtext("Final Model", side=3)
  
  magimage(FinalModel-image, stretchscale = stretchscale, lo = -maximg, 
           hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
           xlab = "x/pix", ylab = "y/pix")
  mtext("Residual", side=3)
  # magimage(ModelNoise, stretchscale = stretchscale, lo = -maximg, 
  #          hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
  #          xlab = "x/pix", ylab = "y/pix")
  # mtext("Final Model+Noise", side=3) 
  
  # magimage(image-FinalModel, stretchscale = stretchscale, lo = -maximg, 
  #          hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
  #          xlab = "x/pix", ylab = "y/pix")
  # mtext("Residual", side=3) 
  
}


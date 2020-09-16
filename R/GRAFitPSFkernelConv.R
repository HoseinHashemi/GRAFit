GRAFitPSFkernelConv <- function(image = image, plot = TRUE) {

  library(magick)
  library(astro)
  library(magicaxis)
  library(OpenImageR)
  # psf_2 = read.fits('~/Desktop/Hosein/HST/Runs/PSF_test/subsample/subsam_psf_.0500.fits')
  # PSFhdr = psf_2$hdr
  
  kern <- matrix(0, ncol = 3, nrow = 3)
  kern[1,1] <- 0.016687; kern[1,2] <- 0.073696; kern[1,3] <- 0.016687
  kern[2,1] <- 0.073696; kern[2,2] <- 0.638469; kern[2,3] <- 0.073696
  kern[3,1] <- 0.016687; kern[3,2] <- 0.073696; kern[3,3] <- 0.016687
  
  conv_PSF <- convolution(image, kernel = kern)
  if (plot) magimage(conv_PSF)
  
  return(conv_PSF)

}

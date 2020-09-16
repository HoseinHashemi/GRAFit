GRAFitPSFMakeStack <- function(psf1, psf2, psf3, psf4) {
  
  error = "ERROR in: PSF downsampling."
  psf1$imDat = profitDownsample(psf1$imDat,3)
  psf2$imDat = profitDownsample(psf2$imDat,3)
  psf3$imDat = profitDownsample(psf3$imDat,3)
  psf4$imDat = profitDownsample(psf4$imDat,3)
  
  # Charge diffusion kernel convolution
  # psf1 = GRAFitPSFkernelConv(psf1, plot = FALSE)
  # psf2 = GRAFitPSFkernelConv(psf2, plot = FALSE)
  # psf3 = GRAFitPSFkernelConv(psf3, plot = FALSE)
  # psf4 = GRAFitPSFkernelConv(psf4, plot = FALSE)
  # error = "ERROR in: PSF charge diffusion kernel convolution."
  # kern <- matrix(0, ncol = 3, nrow = 3)
  # kern[1,1] <- 0.016687; kern[1,2] <- 0.073696; kern[1,3] <- 0.016687
  # kern[2,1] <- 0.073696; kern[2,2] <- 0.638469; kern[2,3] <- 0.073696
  # kern[3,1] <- 0.016687; kern[3,2] <- 0.073696; kern[3,3] <- 0.016687
  
  error = "ERROR in: PSF charge diffusion kernel convolution."
  
  kern1 <- psf1$header[grep(pattern = "COMMENT", psf1$header)]
  kern2 <- psf2$header[grep(pattern = "COMMENT", psf2$header)]
  kern3 <- psf3$header[grep(pattern = "COMMENT", psf3$header)]
  kern4 <- psf4$header[grep(pattern = "COMMENT", psf4$header)]
  
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }
  
  kern1 <- Numextract(kern1)
  kern2 <- Numextract(kern2)
  kern3 <- Numextract(kern3)
  kern4 <- Numextract(kern4)
  kern1 <- matrix(as.numeric(Numextract(kern1)), nrow = 3, ncol = 3)
  kern2 <- matrix(as.numeric(Numextract(kern2)), nrow = 3, ncol = 3)
  kern3 <- matrix(as.numeric(Numextract(kern3)), nrow = 3, ncol = 3)
  kern4 <- matrix(as.numeric(Numextract(kern4)), nrow = 3, ncol = 3)
  
  
  
  psf1 = profitBruteConv(image = psf1$imDat, psf = kern1)
  psf2 = profitBruteConv(image = psf2$imDat, psf = kern2)
  psf3 = profitBruteConv(image = psf3$imDat, psf = kern3)
  psf4 = profitBruteConv(image = psf4$imDat, psf = kern4)
  
  stck_psf = psf1+psf2+psf3+psf4
  finalPSF = stck_psf/sum(stck_psf)
  
  return(list(finalPSF = finalPSF, error = error))
  
}
GRAFitMainFinder = function(src_list = NULL, imDim = NULL, main_pin = FALSE) {
  nSrc = length(src_list$segID) # number of detected sources.
  
  x0 = imDim[1]/2; y0 = imDim[2]/2 # position of cut-out centre.
  sepArr = array(0, dim=nSrc) # empty separation array
  
  # Calculate separation of sources to the centres of image.
  for (i in seq(1,nSrc)) {
    sepArr[i] = sqrt( (src_list$xcen[i] - x0)^2 + (src_list$ycen[i] - y0)^2 )
  }
  
  mainObj = which.min(sepArr) # Main source is the one with the smallest separation from the image centre
  main = src_list[mainObj, ] 
  
  if (main_pin) points(main$xmax, main$ymax, pch= 4, col = 'red')
  
  return( main )
}

# END
  
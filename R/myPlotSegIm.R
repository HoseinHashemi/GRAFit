myPlotSegIm = function(image, segim, mask, sky = 0, main=NULL, ...) 
{
  cmap = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
  # medimg = median(abs(image[image > 0 & Data$region]))
  # maximg = max(abs(image[Data$region]))
  # stretchscale = 1/medimg
  # image = image - sky
  temp = magimage(image, ...)
  if (min(segim, na.rm = TRUE) != 0) {
    segim = segim - min(segim, na.rm = TRUE)
  }
  segvec = which(tabulate(segim) > 0)
  for (i in segvec) {
    z = segim == i
    z = z[ceiling(temp$x), ceiling(temp$y)]
    if (!is.null(main) && i == main){
      contour(temp$x, temp$y, z, add = T, col = rainbow(1000,start=0.0,end=0.01)[sample(1000, 1)], zlim = c(0, 1), lwd=1.5, drawlabels = FALSE, nlevels = 1)
    } else {
      contour(temp$x, temp$y, z, add = T, col = rainbow(1000,start=0.49,end=0.51)[sample(1000, 1)], zlim = c(0, 1), lwd=1., drawlabels = FALSE, nlevels = 1)
    }
    
  }
  if (!missing(mask)) {
    magimage(mask, lo = 0, hi = 1, col = c(NA, hsv(alpha = 0.3)), 
             add = T)
  }
}



### END


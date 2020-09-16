# GRAFitCutout takes the image, ra and dec and the cut rad and returns the cut out image as cut_sci.fits.
# The default cutrad is 400*400 pixel.

#-----------------------
# credit: Hosein Hashemi

GRAFitCutout <- function(dir=input_dir, im=im_sci, Ra=NULL, Dec=NULL, 
                  xrad=400, yrad=400, header= NULL, plot=TRUE){
  
  # xy <- radec2xy(ra_deg, dec_deg, CRVAL1=hdr_info[1], CRPIX1=hdr_info[3], 
  #                CRVAL2=hdr_info[2], CRPIX2=hdr_info[4], CD1_1=hdr_info[5],
  #                CD1_2=hdr_info[6], CD2_1=hdr_info[7], CD2_2=hdr_info[8])
  xy <- magWCSradec2xy(ra_deg, dec_deg, header = header)
  # cut=20/3600       # 20*20 acrsec cut-out. By deviding by 3600 we get degree.
  # sc_pix=c(cut/abs(xscale),cut/abs(yscale))     # the cut out region in pixel.
  # sc_pix=c(xrad,yrad)
  
  x1_0=floor(xy[1] - (xrad/2))
  x2_0=ceiling(xy[1] + (xrad/2))
  y1_0=floor(xy[2] - (yrad/2))
  y2_0=ceiling(xy[2] + (yrad/2))
  
  # read the fits file and cut out the object
  image = read.fits(im, xlo = round(x1_0), xhi = round(x2_0), 
                    ylo = round(y1_0), yhi = round(y2_0))

  if (plot) magimageWCS(image$dat[[1]],image$hdr[[1]],grid.lwd=1 ) 
  
  c_sci=paste(input_dir,'/cut_sci.fits',sep="")
  write.fits(image, file=c_sci)
  
  return(image)
    
}



## COMMENTS ########################
#  c_wht=paste(input_dir,'/cut_wht.fits',sep="")

#  if (im == im_sci) {
#  magimage(image$dat[[1]], xlab= "x [ pixel ]", ylab= "y [ pixel ]")
#write.fits(image, file=c_sci)
#  } else write.fits(image, file=c_wht)

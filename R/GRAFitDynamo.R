# Author: Hosein Hashemi

GRAFitDynamo <- function( im = NULL, loc = c(ra_deg, dec_deg), init_xrad = 300, init_yrad = 300, CutSteps = 100,
                          StarPinPoint = FALSE, blur_sigma = 0.5, GRAFitlib = GRAFitlib, header = NULL, addpad = 0, 
                          pix_scale = NULL, verbose = TRUE, magzero = NULL, tolerance = 5, cut_frc = 0.3, 
                          smooth = TRUE, sigma = 6, sky, pixcut = 5, skycut = 0.2, boundstats = TRUE, 
                          rotstats = TRUE, ImPlot = TRUE, segPlot = FALSE, ... ) {
  
  source(paste(GRAFitlib,'/GRAFitMainFinder.R', sep = ''))
  source(paste(GRAFitlib,'/GRAFitmagcutoutWCS.R', sep=''))
  source(paste(GRAFitlib,'/GRAFitStarPinPoint.R', sep=''))
  plot.new()
  
  xrad = init_xrad; yrad = init_yrad
  main_src_flag_border = 1
  main_src_maj = xrad
  iter = 0      
  
  while( main_src_flag_border != 0 | main_src_maj*4 >= cut_frc*sqrt(2)*xrad ) {  # (main_src$maj*2= kron radius)*2= Kron diameter ---> that we decide Kron diameter of the main source to be smaller than *cut_frac* the cut-out diagonal. to be 
    tryCatch({
      
      xrad = xrad + CutSteps ; yrad = yrad + CutSteps
      iter = iter + 1
      if(verbose) cat( paste( "*** Dynamic cut-out iteration: ", iter, sep = "" ),'\n' )
      # if (iter != 1 ) image=image+sky1
      # cut_im = GRAFitCutout(output_dir, im = im_sci, Ra= ra_deg, Dec= dec_deg,
      #                   xrad = xrad, yrad = yrad, header=sci_hdr, plot= T)
      
      cut_im = GRAFitmagcutoutWCS( image = im, loc = loc, box = c( xrad * pix_scale, yrad * pix_scale ),
                                   paddim = TRUE, addpad = addpad, header = header, grid.lwd = 1, plot = ImPlot )
      
      log3=" Image cutout:: DONE :D "
      image <- cut_im$image
      
      # if (is.na(image) ) main_src_flag_border = 0; break
      
      #----------------------------------------------------------------
      # c_sci=paste(wrk_dir,'/','input_fits','/','cut_sci.fits',sep="")
      
      # segim and sky estimation by ProFit --------------------------
      # seg <- profitMakeSegim(image, tolerance = 10, sigma = 2,
      #                        smooth = TRUE , pixcut = 15, skycut = 2, magzero=ZP,
      #                        pixscale=pix_scale, plot = T,gain = gain_ADU, header = sci_hdr)
      # image= profoundImBlur(image, sigma = 10 , plot = T)
      
      seg <- profoundProFound( image, tolerance = tolerance, sigma = sigma, sky = sky,     # mask = image==0,
                               smooth = smooth , pixcut = pixcut, skycut = skycut, 
                               magzero = magzero, pixscale = pix_scale,
                               header = header, boundstats = boundstats, 
                               rotstats = rotstats, plot = FALSE, ... )  #gain = gain_ADU
      
      segstats <- seg$segstats
      main_src = GRAFitMainFinder(src_list =  segstats, imDim = dim(image))
      main_srcSky = main_src$sky_mean; main_srcSkyRMS = main_src$skyRMS_mean
      
      if (StarPinPoint) {
        
        str_pin_seg = GRAFitStarPinPoint( image = image, mainSeg = seg$segim, blur_sigma = 6, 
                                          segPlot = segPlot, pixscale = pix_scale, 
                                          magzero = magzero, header = header, tolerance = 3, 
                                          sigma = 6, smooth = smooth , pixcut = pixcut, skycut = 1 )
        
        seg$segim = str_pin_seg$finalSeg
      }
      
      # if (iter == 1) sky1 = seg$sky
      
      # image = image - seg$sky
      
      seg_dilate = profoundMakeSegimDilate(image, seg$segim, size = 41, expand = main_src$segID, boundstats = TRUE, 
                                           rotstats = TRUE, header = header, pixscale = pix_scale, magzero = magzero, plot = FALSE )
      
      # if(verbose) cat(" dilation:: DONE :D ",'\n')
      log004=" dilation:: DONE :D "
      segim <- seg_dilate$segim
      mask = segim == 0
      segstats <- seg_dilate$segstats
      main_src = GRAFitMainFinder(src_list =  segstats, imDim = dim(image))
      main_src_flag_border = main_src$flag_border        
      main_src_maj = main_src$semimaj
      
      # loc = c(main_src$RAcen, main_src$Deccen)
    } , error=function(e){
      if (is.na(image) | image[dim(image)/2,dim(image)/2]==0 ) break
      init_xrad = init_xrad+100; init_yrad = init_yrad+100
      # if (all( abs(image - mean(image)) < .001 ) ) break
    }) # END of tryCatch
    
  }  # end of while loop for dynamic cutout.
  
  output = return(list( sky_red_image = image, seg = seg, seg_dilate = seg_dilate, mask = mask, main_src = main_src ))
  
}

# END


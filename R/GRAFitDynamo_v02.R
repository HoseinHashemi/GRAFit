# Author: Hosein Hashemi


GRAFitDynamo_v02 <- function( im = NULL, loc = c(ra_deg, dec_deg), init_xrad = NULL, init_yrad = NULL, 
                             CutSteps = 100, R90 = NULL, StarPinPoint = FALSE, blur_sigma = 0.5, 
                             GRAFitlib = GRAFitlib, header = NULL, pix_scale = NULL, verbose = TRUE, 
                             magzero = NULL, tolerance = 7, cut_frc = 0.3, smooth = TRUE, sigma = 7, 
                             sky, pixcut = 3, skycut = 1.1, boundstats = TRUE, rotstats = TRUE, 
                             ImPlot = TRUE, segPlot = FALSE, output_dir = NULL, ... ) {

  if(verbose) cat( paste( "*** Doing Dynamic cut-out ", sep = "" ),'\n' )
  source(paste(GRAFitlib,'/GRAFitMainFinder.R',sep = ''))
  source(paste(GRAFitlib,'/GRAFitcutoutWCS.R',sep=''))
  source(paste(GRAFitlib,'/GRAFitcutout.R',sep=''))
  source(paste(GRAFitlib,'/GRAFitStarPinPoint.R',sep=''))
  # plot.new()
  
  if (!is.null(init_xrad)) {
    xrad = init_xrad; yrad = init_yrad  # in arcsec  
  } else {
    rad = round(15*R90)   # in arcsec
  }

  xy = magWCSradec2xy(RA = loc[1], Dec = loc[2], header = header)

  cut_im =  GRAFitcutout( GRAFitlib_path = GRAFitlib, image = im, 
                          loc = c(xy[1],xy[2]), box = c( ceiling(rad/0.03), ceiling(rad/0.03) ),
                          paddim = TRUE, header = header, grid.lwd = 1, plot = ImPlot, padVal = 0 )
  
  # cut_im = GRAFitcutoutWCS( GRAFitlib_path = GRAFitlib, image = image, loc = loc, box = c( xrad, yrad ),
  #                              paddim = TRUE, header = header, grid.lwd = 1, plot = ImPlot, padVal = 0 )
  # 
  
  image <- cut_im$image

  invisible(
  ifelse ( (image[1,] == 0) | ( image[,1] == 0 ), 
           image <- floodFill( cut_im$image, c(1,1), col = NA ), 
           (image = image) )
  )
  invisible(
  ifelse ( ( image[dim(image)[1], ] == 0 ) | (image[ ,dim(image)[2] ] == 0 ), 
           image <- floodFill( cut_im$image, c(dim(image)[1], dim(image)[2]), col = NA ), 
           (image = image) ) 
  )
  
  # mask = matrix(0, dim(image)[1], dim(image)[1]); mask[which(is.na(image))] = 1
  # mask = matrix(0, dim(image)[1], dim(image)[1]); mask[which(image == 0)] = 1
  box_car = round(5*R90/pix_scale)
  
  print("----- profound") 
  print(paste("box = ", box_car, "R90 = ", R90))

  seg <- profoundProFound( image = image, tolerance = tolerance, sigma = sigma, sky = sky,      # mask = image==0,
                           smooth = smooth , pixcut = pixcut, skycut = skycut,
                           magzero = magzero, pixscale = pix_scale, 
                           box = c(box_car,box_car), type = "bicubic", 
                           header = header, boundstats = boundstats, groupstats = FALSE, 
                           rotstats = rotstats, plot = FALSE, ... )  #gain = gain_ADU

  segstats <- seg$segstats
  main_src0 = GRAFitMainFinder(src_list =  segstats, imDim = dim(image))
  
  print("----- dilation")  
  
  seg_dilate = profoundMakeSegimDilate(image = image, seg$segim, size = 51, #expand = main_src$segID,
                                       boundstats = TRUE, rotstats = TRUE, header = header,
                                       pixscale = pix_scale, magzero = magzero, plot = FALSE )

# seg_dilate = profoundMakeSegimExpand(image = image, seg$segim, expand = main_src$segID,
#                                      skycut = -1,
#                                      boundstats = TRUE, rotstats = TRUE, header = header, 
#                                      pixscale = pix_scale, magzero = magzero, plot = FALSE )

  log004 = " dilation:: DONE :D "
  segim_dil <- seg_dilate$segim
  mask = segim_dil == 0
  segstats <- seg_dilate$segstats
  main_src = GRAFitMainFinder(src_list =  segstats, imDim = dim(image))

  i = 5
  # flag_border = 1
  # while( flag_border != 0 ) {
  #   i = i + 1
    box = c(round(i*R90/pix_scale), round(i*R90/pix_scale))

    # cut_im_final_skyred <- GRAFitcutout(GRAFitlib_path = GRAFitlib, image = seg$image, 
    #                                     loc = c(main_src$xmax, main_src$ymax), box = box, plot = ImPlot)
    # cut_im_final_orig <- GRAFitcutout( GRAFitlib_path = GRAFitlib, image = image,
    #                                   loc = c(main_src$xmax, main_src$ymax), box = box )
    cut_im_final_orig <- magcutoutWCS( image = im, header = header, 
                                       loc = loc, 
                                       box = c(round(i*R90), round(i*R90)) )

    # loc = cut_im$loc
    loc = c(main_src$xcen, main_src$ycen)
    box = dim(cut_im_final_orig$image)
      
    segim <- GRAFitcutout( GRAFitlib_path = GRAFitlib, image = segim_dil, 
                              box = box )
    
    sky <- magcutout( seg$sky, box = box )
    
    skyRMS <- magcutout( seg$skyRMS, box = box )
    
    objects <- magcutout( seg_dilate$objects, box = box )
    
    box_car = round(3*main_src$R90/pix_scale)
    seg2 <- profoundProFound( image = cut_im_final_orig$image, segim = segim$image , 
                              tolerance = tolerance, sigma = sigma,     # mask = image==0,
                              smooth = smooth , pixcut = pixcut, skycut = skycut,
                              magzero = magzero, pixscale = pix_scale, box = c(box_car, box_car),
                              header = header, boundstats = boundstats, groupstats = FALSE,
                              rotstats = rotstats, plot = ImPlot, lowmemory = F, ... )
    
    main_src2 = GRAFitMainFinder(src_list =  seg2$segstats, imDim = dim(cut_im_final_orig$image))
    # flag_border <- main_src2$flag_border
  # }
  
  if (!is.null(output_dir)) {
    CairoPNG(filename = paste(output_dir,"/postage_stamps.png", sep = ''), 
             width = 1000, height = 1068)
      par(mfrow=c(2,2), cex.lab = 2, cex.axis =2)
      profoundSegimPlot(image, seg$segim)
      title(main = "Initial cut out = 15R90 (Y UltraVISTA)", cex.main = 2)
      profoundSegimPlot(cut_im_final_orig$image, segim$image)
      title(main = paste("Final cut out = ",i, "R90 (Y UltraVISTA)", sep = ''), cex.main = 2)
      magimage(seg$sky)
      title(main = "sky", cex.main = 2)
      profoundSkyEst(image = image, objects = seg$objects_redo, plot = TRUE)
      title(main = "sky pix val distribution", cex.main = 2)
    dev.off()
    
    log3 = " Image cutout:: DONE :D "
    output = return(list( orig_image = cut_im_final_orig$image, sky_red_image = cut_im_final_orig$image-sky$image, 
                          segim = segim$image, main_src = main_src2, sky = sky$image, skyRMS = skyRMS$image,
                          objects = objects$image))
    
  } else {
    log3 = " Image cutout:: DONE :D "
    output = return(list( orig_image = cut_im_final_orig$image, sky_red_image = cut_im_final_orig$image-sky$image, 
                          segim = segim$image, main_src = main_src2, sky = sky$image, skyRMS = skyRMS$image,
                          objects = objects$image))
    
  }
  
}

# END


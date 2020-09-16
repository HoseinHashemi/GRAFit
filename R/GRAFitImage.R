GRAFitImage <- function(ID = NULL, loc = c(RA, DEC), plot = TRUE, AddSeg = TRUE,
                        data_dir = data_dir, GRAFitlib = GRAFitlib, ...) {

  source(paste(GRAFitlib,'/GRAFitDynamo_v2.R',sep=''))
  source(paste(GRAFitlib,'/GRAFitcutoutWCS.R',sep=''))
  source(paste(GRAFitlib,'/GRAFitcutout.R',sep=''))
  source(paste(GRAFitlib,'/GRAFitMainFinder.R',sep=''))
  source(paste(GRAFitlib,'/GRAFitFrameFinder.R',sep=''))
  library(data.table)  
  
  D10 = suppressWarnings(fread(file = paste(data_dir,'/D10_XY_shuf.csv', sep = '')))
  if (!is.null(ID)) {
    ra_deg = D10[D10$D10CATAID == ID, RA]
    dec_deg = D10[D10$D10CATAID == ID, DEC]
    R90 = D10[D10$D10CATAID == ID, R90]
    redshift = D10[D10$D10CATAID == ID, ZBEST]
    stellar_mass = D10[D10$D10CATAID == ID, STELLARMASS ]
  } else {
    ra_deg = loc[1]
    dec_deg = loc[2]
  }
  
  ra_form = as.numeric(format(ra_deg, digits = 8))
  dec_form = as.numeric(format(dec_deg, digits = 8))
  
  frame_im_name = GRAFitFrameFinder(GRAFitlib = GRAFitlib, data_dir = data_dir,
                                    target_loc = c(ra_deg, dec_deg))
  
  frame_im = readFITS( frame_im_name )
  frame_hdr = frame_im$hdr
  frame_im= frame_im$imDat
  
  # Read required info from header.
  PHOTFLAM = as.numeric( frame_hdr[which(frame_hdr == "PHOTFLAM")+1] )
  PHOTPLAM = as.numeric( frame_hdr[which(frame_hdr == "PHOTPLAM")+1] )
  PHOTZPT = as.numeric( frame_hdr[which(frame_hdr == "PHOTZPT")+1] )
  EXPTIME = as.numeric( frame_hdr[which(frame_hdr == "EXPTIME")+1] )
  PHOTBW = as.numeric( frame_hdr[which(frame_hdr == "PHOTBW")+1] )
  CCDGAIN = as.numeric( frame_hdr[which(frame_hdr == "CCDGAIN")+1] )
  FILTER = as.character( frame_hdr[which(frame_hdr == "FILTER2")+1] )
  # pix_scale = profoundGetPixScale( frame_hdr )
  pix_scale = 0.03

  gain_ADU <- EXPTIME * CCDGAIN    # gain_ADU = EXPTIME * CCD_GAIN (=1)
  frame_elec_map <- gain_ADU * frame_im

  # ZERO POINT ####################
  instr_ABMAG_ZPT = -2.5 * log10(PHOTFLAM) - 5*log10(PHOTPLAM)-2.408
  ZP = instr_ABMAG_ZPT + 2.5*log10(gain_ADU)
  
  #----Dynamic cut-out----------------------------------------------------
  error = "ERROR in: Dynamoic cut out"
  dyn_cut = GRAFitDynamo_v2( image = frame_elec_map, loc = c(ra_deg, dec_deg), R90 = R90,
                             GRAFitlib = GRAFitlib, gain = 1, 
                             magzero = ZP, pix_scale = pix_scale, header = frame_hdr, 
                             ImPlot = FALSE, segPlot = FALSE, output_dir = output_dir )
  
  main_src = dyn_cut$main_src
  image = dyn_cut$sky_red_image
  segim = dyn_cut$segim
  mask = segim == 0
  
  if (plot) {
    cmap = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
    medimg = median(abs(image[image > 0 & dyn_cut$segim == main_src$segID]))
    maximg = max(abs(image[dyn_cut$segim == main_src$segID]))
    stretchscale = 1/medimg

      image = magimage(image, stretchscale = stretchscale, lo = -maximg,
                       hi = maximg, type = "num", zlim = c(0, 1), col = cmap, ...)  

    
    if ( AddSeg ) {
      segvec = which(tabulate(segim) > 0)
      col = rainbow(max(segim), 
                    end = 2/3)
      for (i in segvec) {
        z = segim == i
        z = z[ceiling(image$x), ceiling(image$y)]
        contour(image$x, image$y, z, add = T, col = col[i], lwd = 2,
                zlim = c(0, 1), drawlabels = FALSE, nlevels = 1)
      }  
    }
    
    # image = profoundSegimPlot(image = image, segim = segim, stretchscale = stretchscale, lo = -maximg,
    #                                    hi = maximg, type = "num", zlim = c(0, 1), col = cmap, ... )
    
  }

  return(list(SkyRedImage = image, redshift = redshift, stellar_mass = stellar_mass))
  
}

  
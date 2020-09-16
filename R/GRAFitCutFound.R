GRAFitCutFound <- function(wrk_dir=NULL, data_dir= NULL, object_list=NULL, ncores=1, verbose=TRUE ) {
  
  # .libPaths(c("/home/arobotham/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))
  options(digits=10)
  assign("wrk_dir", wrk_dir, envir = .GlobalEnv) 
  assign("data_dir", data_dir, envir = .GlobalEnv) 
  
  # wrk_dir='/gamah/hosein/HST/100_sample3/MCMC'
  # data_dir='/gamah/hosein/HST/data/COSMOS'
  # R_func_dir='/gamah/hosein/HST/R_func'
  # catalog_dir='/gamah/hosein/HST/data/catalogs'
  #   R_func_dir='/Users/22111305/Desktop/PhD/HST/code/R_func'
  #   catalog_dir='/Users/22111305/Desktop/PhD/HST/data/catalogs'
  #Robin's laptop
  R_func_dir= '/Users/robincook/Desktop/Hosein/HST/code/R_func'
  catalog_dir='/Users/robincook/Desktop/Hosein/HST/data/catalogs'
  # Load the required libraries. 
  library(knitr) 
  library(ProFit)
  library(FITSio)
  library(ProFound)
  library(astro)
  library(celestial)
  library(devtools)
  library(magicaxis)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(LaplacesDemon)
  library(RColorBrewer)
  registerDoParallel(cores = ncores)     # snow-like funtionality.
  
  source(paste(R_func_dir,'/HST_cut.r',sep=''))
  source(paste(R_func_dir,'/HST_ap.r',sep=''))
  source(paste(R_func_dir,'/TinyTimACSR.R',sep=''))
  source(paste(R_func_dir,'/TinyTimACSR_labtop.R',sep=''))
  source(paste(R_func_dir,'/TinyTimACSR_Robin_lab.R',sep=''))
  source(paste(R_func_dir,'/HST_mu.r',sep=''))
  source(paste(R_func_dir,'/FITS_hdr_list.R',sep=''))
  source(paste(R_func_dir,'/add_pseudo_bulge.R',sep=''))
  source(paste(R_func_dir,'/myPlotSegIm.R',sep=''))
  source(paste(R_func_dir,'/mainFinder.R',sep=''))
  
  # set working directory --------------
  setwd(wrk_dir)
  
  # ######################################  
  # ########## READ DATA FITS HEADERS ####
  # ######################################
  # Read FITS files and write the TARGET'S RA AND DEC (i.e. the centre of frames)
  # in a file name: target_cor.csv  ----------------------------
  #  FITS_hdr_list(data_path = data_dir)   # This function reads all the fits header and save them in a file named: target_cor.csv
  file_names = list.files(path = data_dir ,full.names = TRUE, pattern="*.fits")
  
  trgt_cor = read.csv(file = paste(data_dir,'/target_cor.csv',sep = ""),sep = ",")
  RA_TARG <- trgt_cor[,1]
  DEC_TARG <- trgt_cor[,2]
  trgt_cor<-0
  
  out_catalog= paste(wrk_dir,'/','catalog.csv',sep = "")
  CATAID <- object_list$CATAID
  RA <- object_list$RA
  DEC <- object_list$DEC
  z <- object_list$Z_BEST
  LBT <- object_list$look_bk_T
  maj_ax <- object_list$SEMIMAJ_AS
  Ymag <- object_list$Ymag
  stellar_mass <- object_list$mass_stellar_best_fit
  
  ######################################  
  ########## PREPARE INPUT FOR PROFIT ##########
  ######################################
  system(paste('mkdir ', wrk_dir,'/cut_and_seg500test', sep=''))
  output_dir = paste(wrk_dir,'/cut_and_seg500test', sep='')
  print("step1")
  foreach (j=1:NROW(object_list)) %dopar% {     # Loop on all objects. ,.errorhandling="pass" , .combine = 'cbind', .multicombine = TRUE
    
    tryCatch({
      #       source(paste(R_func_dir,'/TinyTimACSR_labtop.R',sep=''))
      source(paste(R_func_dir,'/TinyTimACSR.R',sep=''))     
      source(paste(R_func_dir,'/add_pseudo_bulge.R',sep=''))
      source(paste(R_func_dir,'/HST_mu.r',sep=''))
      source(paste(R_func_dir,'/TinyTimACSR_Robin_lab.R',sep=''))
      ###############################################################
      
      t1 <- proc.time()       
      ra_deg=RA[j]
      dec_deg=DEC[j]
      print(c(ra_deg,dec_deg))
      log01=paste(' CATAID :: ',CATAID[j], sep = '')
      log1=paste(' RA & DEC :: ',ra_deg,'+',dec_deg, sep = '') 
      log001=paste(' Redshift:: ', z[j])
      # ra_form=as.numeric(format(ra_deg, digits = 8))
      # dec_form=as.numeric(format(dec_deg, digits = 8))
      
      CATAID_out=CATAID[j];z_out=z[j];look_back_t_out=LBT[j];RA_out=ra_deg;DEC_out=dec_deg
      LAMBDAR_maj_ax_out=maj_ax[j];Ymag_out=Ymag[j];stellar_mass_out=stellar_mass[j]
      
      ang_sep <- vector()
      foreach (i=1:length(file_names)) %do% {
        ang_sep[i] = acos( sin(dec_deg)*sin(DEC_TARG[i])+cos(dec_deg)*cos(DEC_TARG[i])*cos(ra_deg-RA_TARG[i]) )
      }
      
      im_sci = file_names[which.min(ang_sep)]
      if(verbose)cat(' Found galaxy in frame:: ',im_sci, '\n')
      log2=paste(' Found galaxy in frame:: ',im_sci)
      
      sci=readFITS(im_sci)
      sci_hdr = sci$hdr
      im= sci$imDat
      
      PHOTFLAM = as.numeric(sci_hdr[which(sci_hdr=="PHOTFLAM")+1])
      PHOTPLAM = as.numeric(sci_hdr[which(sci_hdr=="PHOTPLAM")+1])
      PHOTZPT = as.numeric(sci_hdr[which(sci_hdr=="PHOTZPT")+1])
      EXPTIME = as.numeric(sci_hdr[which(sci_hdr=="EXPTIME")+1])
      PHOTBW = as.numeric(sci_hdr[which(sci_hdr=="PHOTBW")+1])
      #    CCDGAIN = as.numeric(sci_hdr[which(sci_hdr=="CCDGAIN")+1]) 
      CCDGAIN=1
      EXPSTART = as.numeric(sci_hdr[which(sci_hdr=="EXPSTART")+1])
      CCDCHIP= as.numeric(sci_hdr[which(sci_hdr=="CCDCHIP")+1])
      FILTER = as.character(sci_hdr[which(sci_hdr=="FILTER2")+1])
      DATE_OBS = as.character(sci_hdr[which(sci_hdr=="DATE-OBS")+1])
      PIXELSCALE = 0.05
      
      #-------------------------------------------------
      # 400*400 pixel cutout.---------------------------
      xrad = 400; yrad = 400
      #profitGetPixScale(sci_hdr)
      # cut_im = HST_cut(output_dir, im = im_sci, Ra= ra_deg, Dec= dec_deg,
      #                   xrad = xrad, yrad = yrad, header=sci_hdr, plot= T)
      
      cut_im=magcutoutWCS(image = im, loc = c(ra_deg, dec_deg), box = c(xrad*0.03,yrad*0.03),
                          header = sci_hdr, grid.lwd=1, plot= FALSE)
      
      # write.fits(cut_im$image, file=paste(output_dir,'/cut_im.fits',sep=""))
      # cut_hdr <- read.fitshdr(im_sci)
      # write.fitshdr(cut_hdr, file=paste(output_dir,'/cut_im.fits',sep=""))

      log3=" Image cutout:: DONE :D "
      
      #----------------------------------------------------------------
      # c_sci=paste(wrk_dir,'/','input_fits','/','cut_sci.fits',sep="")
      
      # make thes electron map ############################
      gain_ADU <- EXPTIME * CCDGAIN    # gain_ADU = EXPTIME * CCD_GAIN (=1)
      electron_map <- gain_ADU * cut_im$image
      #elec_mp=paste(wrk_dir,'/','input_fits','/','electron_map.fits',sep="")
      #write.fits(electron_map, file=elec_mp)
      image <- electron_map
      # ZERO POINT ####################
      instr_ABMAG_ZPT = -2.5 * log10(PHOTFLAM) - 5*log10(PHOTPLAM)-2.408
      #ZP = instr_ABMAG_ZPT + 2.5*log10(multiplier)
      ZP = instr_ABMAG_ZPT + 2.5*log10(gain_ADU)
      
      #----------------------------------------------------------------
      print("step3")
      # segim and skyRMS by ProFit --------------------------
      # seg <- profitMakeSegim(image, tolerance = 10, sigma = 2,
      #                        smooth = TRUE , pixcut = 15, skycut = 2, magzero=ZP,
      #                        pixscale=PIXELSCALE, plot = T,gain = gain_ADU, header = sci_hdr)
      
      seg <- profitProFound(image, tolerance = 5, sigma = 3, # mask = image==0,
                            smooth = TRUE , pixcut = 5, skycut =2, magzero=ZP,  
                            gain = gain_ADU, header = sci_hdr, boundstats = TRUE, 
                            rotstats = TRUE , size = 11, plot = T)
      
      log4=" Segmentation Map:: DONE :D "
      image = image - seg$sky
      segstats <- seg$segstats
      
      main_src=mainFinder(segstats,dims = dim(image))
      
      if (NROW(main_src) != 0) {
        log04=" Main source finding:: DONE :D "
      } else {
        log04= " No main source found "
        break
      }
      
      seg_dilate = profoundMakeSegimDilate(image, seg$segim, plot=T, size= 85,expand=main_src$segID,
                                           header = sci_hdr,pixscale = PIXELSCALE, magzero = ZP)
      if(verbose)cat("dilation:: DONE :D ")
      log004="dilation:: DONE :D "
      segim <- seg_dilate$segim
      mask=segim==0
      segstats <- seg_dilate$segstats
      print("step4")
      profitSegimPlot(image, segim, sky = seg$sky, main=main_src$segID,xlab="x/pix",ylab="y/pix",mask = mask)
      # save the segmentation map as a png file. ---------------------
      png(file=paste(output_dir,"/G",CATAID[j],".png",sep=''),width=10.8,height=4,units="in",res=200) #,width=7,height=7.6,units="in",res=150)
      par(mfrow=c(1,3))
      magimage(cut_im$image,xlab="x/pix",ylab="y/pix")
      mtext("Image", side=3, line = 1)
      myPlotSegIm(image, segim, sky = seg$sky, main=main_src$segID,xlab="x/pix")
      mtext("Sky-subtracted Image+sources", side=3, line = 1)
      magimage(segim, col=c(0,rainbow(100)),xlab="x/pix")
      mtext("Segmentation Map", side=3, line = 1)
      dev.off()
       # if(verbose)cat('Seg Map saved:: DONE :D ', '\n', file=paste(wrk_dir,'/',logfile, sep = '' ), append=TRUE)
      
      # SIGMA by ProFit ------------------------------------------
      
      sigma=profitMakeSigma(image, objects=seg_dilate$objects, sky=seg$sky,
                            skyRMS =seg$skyRMS, image_units = 'elec', sky_units = 'elec',
                            read_units = 'elec', dark_units = 'elec', output_units = 'elec', plot=T, xlab="x/pix",ylab="y/pix")
magimage(seg$skyRMS, xlab="x/pix",ylab="y/pix")
      # png(file=paste(output_dir,"/G",CATAID[j],"_sigma.png",sep=''),width=5,height=5.6,units="in",res=200)
      # magimage(sigma,xlab="x/pix",ylab="y/pix")
      # dev.off()
      
      log5= " Sigma map:: DONE :D "
      
      
    }, error=function(e){
      print("error")
    })    # End of tryCatch
  }       # End of the loop on galaxies 
  
}

#END
      
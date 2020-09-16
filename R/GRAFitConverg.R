#devtools::install_github("renkun-ken/formattable")
#install_github('asgr/magicaxis',force=T)
#install_github("ICRAR/ProFit")
#install_github('asgr/ProFound')
# source("http://bioconductor.org/biocLite.R")
# biocLite("EBImage")


GRAFitConverg <- function(wrk_dir=NULL, data_dir= NULL, object_list=NULL, ncores=1,
                       logfile="logfile.txt", nComp=2,verbose=TRUE ) {
  if(NROW(object_list) > 1) stop("Only one object per run is accepted.")
  # .libPaths(c("/home/arobotham/R/x86_64-pc-linux-gnu-library/3.2",.libPaths()))
  options(digits=10)
  assign("wrk_dir", wrk_dir, envir = .GlobalEnv) 
  assign("data_dir", data_dir, envir = .GlobalEnv) 
  
  # wrk_dir='/gamah/hosein/HST/100_sample3/MCMC'
  # data_dir='/gamah/hosein/HST/data/COSMOS'
  R_func_dir='/gamah/hosein/HST/R_func'
  catalog_dir='/gamah/hosein/HST/data/catalogs'
#     R_func_dir='/Users/22111305/Desktop/PhD/HST/code/R_func'
#     catalog_dir='/Users/22111305/Desktop/PhD/HST/data/catalogs'
  #Robin's laptop
#   R_func_dir= '/Users/robincook/Desktop/Hosein/HST/code/R_func'
#   catalog_dir='/Users/robincook/Desktop/Hosein/HST/data/catalogs'
  # Load the required libraries. 
  library(knitr) 
  library(ProFit)
  library(FITSio)
  library(ProFound)
  library(astro)
  library(celestial)
  # library(LAMBDAR)
  library(devtools)
  library(magicaxis)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(LaplacesDemon)
  library(RColorBrewer)
  registerDoParallel(cores = ncores)     # snow-like funtionality.
  # getDoParWorkers()
  # parallel::detectCores()
  # cl <- makeCluster(2)
  # cl <- makePSOCKcluster(4)
  # cl <- makeMPIcluster(4)
  # cl <- makeCluster(ncores)
  # library(doSNOW)
  # cl <- makeCluster(ncores)
  # registerDoSNOW(cl)
  # library(doSNOW)
  # cluster = makeCluster(4, type = "MPI")
  # registerDoSNOW(cluster)
  
  source(paste(R_func_dir,'/HST_cut.r',sep=''))
  source(paste(R_func_dir,'/HST_ap.r',sep=''))
  source(paste(R_func_dir,'/TinyTimACSR.R',sep=''))
  source(paste(R_func_dir,'/TinyTimACSR_labtop.R',sep=''))
  source(paste(R_func_dir,'/TinyTimACSR_Robin_lab.R',sep=''))
  source(paste(R_func_dir,'/HST_mu.r',sep=''))
  source(paste(R_func_dir,'/FITS_hdr_list.R',sep=''))
  source(paste(R_func_dir,'/add_pseudo_bulge.R',sep=''))
  source(paste(R_func_dir,'/myPlotSegIm.R',sep=''))
  
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
  # Bulge/Total flux & Bulge effective radius=bulgeRe*ProFoundRe, bulgeRe=c(0.3,0.6,1.3,1.6)        
  BTratio_BRe=list(c(0.3,0.1),c(0.3,0.2),c(0.3,0.3),c(0.3,0.4),c(0.3,0.5),c(0.3,0.6),
               c(0.6,0.1),c(0.6,0.2),c(0.6,0.3),c(0.6,0.4),c(0.6,0.5),c(0.6,0.6),
               c(1.3,0.1),c(1.3,0.2),c(1.3,0.3),c(1.3,0.4),c(1.3,0.5),c(1.3,0.6),
               c(1.6,0.1),c(1.6,0.2),c(1.6,0.3),c(1.6,0.4),c(1.6,0.5),c(1.6,0.6)) 
  
  # foreach (j=1:NROW(object_list)) %do% {     # Loop on all objects. ,.errorhandling="pass" , .combine = 'cbind', .multicombine = TRUE
    j=1
    tryCatch({
      #       source(paste(R_func_dir,'/TinyTimACSR_labtop.R',sep=''))
      source(paste(R_func_dir,'/TinyTimACSR.R',sep=''))     
      source(paste(R_func_dir,'/TinyTimACSR_Robin_lab.R',sep=''))
      ###############################################################
      log01=' CATAID :: NA :( '; log1=' RA & DEC:: NA :( '; log001= 'Redshift:: NA :( '; log2=' Found galaxy in frame:: NA :( '
      log3=' Image cutout:: NA :( '; log4=' Segmentation Map:: NA :( '; log04= " Main source found:: NA ";
      log004="segmentation dilated:: NA :( "; log5=' Sigma map:: NA :( '
      log6=' PSF:: NA :( '; log7=' Initial Model:: NA :( '; log8=' optimfit:: NA :( ';log9=' Laplace Approximation:: NA :( ';
      log09=' Testing LA conversion:: NA :( '; log009='LAfit converged/diverged ? NA :(';
      log11=' Laplaces Demon:: NA :( '; log12=' Profile Plot:: NA :( '; log13=' Output Catalogue:: NA :( '; log14= ' Elapsed time:: NA :( '
      
      CATAID_out=999.9;z_out=999.9;look_back_t_out=999.9;RA_out=999.9;DEC_out=999.9;LAMBDAR_maj_ax_out=999.9;Ymag_out=999.9;stellar_mass_out=999.9
      ProFound_R50_out=999.9;ProFound_R100_out=999.9;ProFound_maj_ax_out=999.9;ProFound_mag_out=999.9;edge_frac_out=999.9;asym_out=999.9
      Init_mag1=999.9;Init_mag2=999.9;Init_re1=999.9
      Init_re2=999.9;Init_n1=999.9;Init_n2=999.9;Init_ang1=999.9;Init_ang2=999.9
      Init_axrat1=999.9;Init_axrat2=999.9
      MCMC_xcen1=999.9;MCMC_xcen2=999.9
      MCMC_ycen1=999.9;MCMC_ycen2=999.9;MCMC_mag1=999.9;MCMC_mag2=999.9;MCMC_re1=999.9
      MCMC_re2=999.9;MCMC_n1=999.9;MCMC_n2=999.9;MCMC_ang1=999.9;MCMC_ang2=999.9
      MCMC_axrat1=999.9;MCMC_axrat2=999.9;MCMC_box1=999.9;MCMC_box2=999.9;elapsed_time=999.9
      ###############################################################
      
      t1 <- proc.time()       
      ra_deg=RA[j]
      dec_deg=DEC[j]
      print(c(ra_deg,dec_deg))
      log01=paste(' CATAID :: ',CATAID[j], sep = '')
      log1=paste(' RA & DEC :: ',ra_deg,'+',dec_deg, sep = '') 
      log001=paste(' Redshift:: ', z[j])
      ra_form=as.numeric(format(ra_deg, digits = 8))
      dec_form=as.numeric(format(dec_deg, digits = 8))
      system(paste('mkdir ', wrk_dir,'/',ra_form,'+',dec_form, sep=''))
      output_dir = paste(wrk_dir,'/',ra_form,'+',dec_form, sep='')
      
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
      # DATE_OBS = as.character(sci_hdr[which(sci_hdr=="DATE-OBS")+1])
      PIXELSCALE = profitGetPixScale(sci_hdr)
      
      #-------------------------------------------------
      # 400*400 pixel cutout.---------------------------
      xrad = 400; yrad = 400
      
      # cut_im = HST_cut(output_dir, im = im_sci, Ra= ra_deg, Dec= dec_deg,
      #                   xrad = xrad, yrad = yrad, header=sci_hdr, plot= T)
      
      cut_im=magcutoutWCS(image = im, loc = c(ra_deg, dec_deg), box = c(xrad*0.03,yrad*0.03),
                          header = sci_hdr, grid.lwd=1, plot= F)
      
      write.fits(cut_im$image, file=paste(output_dir,'/cut_im.fits',sep=""))
      cut_hdr <- read.fitshdr(im_sci)
      write.fitshdr(cut_hdr, file=paste(output_dir,'/cut_im.fits',sep=""))
      
      log3=" Image cutout:: DONE :D "
      
      #----------------------------------------------------------------
      # c_sci=paste(wrk_dir,'/','input_fits','/','cut_sci.fits',sep="")
      
      # make the electron map ############################
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
      
      # segim and skyRMS by ProFit --------------------------
      # seg <- profitMakeSegim(image, tolerance = 10, sigma = 2,
      #                        smooth = TRUE , pixcut = 15, skycut = 2, magzero=ZP,
      #                        pixscale=PIXELSCALE, plot = T,gain = gain_ADU, header = sci_hdr)
      
      seg <- profitProFound(image, tolerance = 5, sigma = 3, # mask = image==0,
                            smooth = TRUE , pixcut = 5, skycut =2, magzero=ZP,  
                            gain = gain_ADU, header = sci_hdr, boundstats = TRUE, 
                            rotstats = TRUE , size = 11, plot = F)
      
      log4=" Segmentation Map:: DONE :D "
      
      segstats <- seg$segstats
      image = image - seg$sky
      main_src = segstats[which(xrad/2-20 < xrad-segstats$xcen & xrad-segstats$xcen < xrad/2+20 &
                                  yrad/2-20 < yrad-segstats$ycen & yrad-segstats$ycen < yrad/2+20 ), ]
      
      if (NROW(main_src) != 0) {
        log04=" Main source finding:: DONE :D "
      } else {
        log04= " No main source found "
        break
      }
      
      seg_dilate = profitMakeSegimDilate(image, seg$segim, plot=F, size= 85,expand=main_src$segID)
      if(verbose)cat("dilation:: DONE :D ")
      log004="dilation:: DONE :D "
      segim <- seg_dilate$segim
      mask=segim==0
      segstats <- seg_dilate$segstats
      
      ProFound_R50_out=main_src$R50;ProFound_R100_out=main_src$R100
      ProFound_maj_ax_out=main_src$maj;ProFound_mag_out=main_src$mag;
      edge_frac_out=main_src$edge_frac;asym_out=main_src$asymm
      
      
      #       myPlotSegIm(image, segim, sky = seg$sky, main=main_src$segID )
      
      # save the segmentation map as a png file. ---------------------
      #       jpeg(filename = paste(output_dir,'/seg.jpeg',sep = '')) #,width=7,height=7.6,units="in",res=150)
      #       c<-profitSegimPlot(image,seg$segim)
      #       dev.off()
      # if(verbose)cat('Seg Map saved:: DONE :D ', '\n', file=paste(wrk_dir,'/',logfile, sep = '' ), append=TRUE)
      
      # SIGMA by ProFit ------------------------------------------
      
      sigma=profitMakeSigma(image, objects=seg_dilate$objects, sky=seg$sky, 
                            skyRMS =seg$skyRMS, image_units = 'elec', sky_units = 'elec',
                            read_units = 'elec', dark_units = 'elec', output_units = 'elec', plot=F)
      
      log5= " Sigma map:: DONE :D "
      
      # Write sigma.fits -----------------------------------------
      # sig=paste(wrk_dir,'/','input_fits','/','sigma.fits',sep="")
      # write.fits(sigma, file=sig)
      
      # estimate PSF using TinyTim -------------------------------
      # find the location of the object on CCD in pixel ---------
      xy <- magWCSradec2xy(ra_deg, dec_deg, header = sci_hdr)
      
      # read the focus values from model. Provided by STScI. 
      focus_tab <- read.table(file = paste(catalog_dir,'/Focus.txt', sep = ""), sep = "", 
                              col.names = c("Julian date","month","day","year","time","focus_val"))
      
      focus <- focus_tab[which.min(abs(EXPSTART-focus_tab$Julian.date)),6]
      
      #     The x and y that are converted to pixel by "radec2xy" 
      #     should be scaled by the ACS pixel scale which is= 0.05
      
      if (object_list$EBV[j] > 0) ebmv=object_list$EBV[j] else ebmv=0   # E(B-V) extinction values calculated from Schlegel+98 dust maps and taken from Laigle+2016
      if(verbose)cat("Making PSF ..." ,'\n')
      # TinyTimACSR_Robin_lab(output_dir, 'psf', CCDCHIP , xy[1]*PIXELSCALE, xy[2]*PIXELSCALE, FILTER, focus, ebmv=ebmv)
      TinyTimACSR(output_dir, 'psf', CCDCHIP , xy[1]*PIXELSCALE, xy[2]*PIXELSCALE, FILTER, focus, ebmv=ebmv)
      
      psf = readFITS(paste(output_dir,'/psf_image00.fits',sep=""))$imDat
      log6=" PSF :: DONE :D "
      
      # ProFit -------------------------------
      # if (nComp == 2) {
      foreach(ii=1:length(BTratio_BRe)) %dopar% {
          ang <- c(0,main_src$ang)
          xcen <- main_src$xcen
          ycen <- main_src$ycen
          # mag <- -2.5*log10(c(main_src$flux*0.1,main_src$flux*0.9))-48.6
          # mag <- main_src$mag
          flux <- main_src$flux
          mag <- -2.5*log10(c(main_src$flux*BTratio_BRe[[ii]][2],main_src$flux*(1-BTratio_BRe[[ii]][2])))+ZP
          re=sqrt(main_src$N50/(pi*main_src$axrat))*BTratio_BRe[[ii]][1]
          axrat<-c(1,(main_src$axrat))
          
          modellist=list(
            sersic=list(
              xcen= rep(xcen,2),
              ycen= rep(ycen,2),
              mag= mag ,
              re= c(0.2*re,re),
              nser= c(4, 1.0000),
              ang= ang, 
              axrat= axrat, #min/maj: 1= o, 0= |
              box=c(0, 0)     
            )
          )
          
          tofit=list(
            sersic=list(
              xcen= c(TRUE,NA), # fit for xcen of both disk and bulge.
              ycen= c(TRUE,NA), # fit for ycen of both disk and bulge.
              mag= c(TRUE,TRUE), # Fit for both
              re= c(TRUE,TRUE), # Fit for both
              nser= c(TRUE,TRUE), # Fit for bulge
              ang= c(FALSE,TRUE), # Fit for disk
              axrat= c(FALSE,TRUE), # Fit for disk
              box= c(FALSE,FALSE) # Fit for neither
            )
          )
          
          # What parameters should be fitted in log space:
          
          tolog=list(
            sersic=list(
              xcen= c(FALSE,FALSE),
              ycen= c(FALSE,FALSE),
              mag= c(FALSE,FALSE),
              re= c(TRUE,TRUE),    #re is best fit in log space
              nser= c(TRUE,TRUE),  #nser is best fit in log space
              ang= c(FALSE,FALSE),
              axrat= c(TRUE,TRUE), #axrat is best fit in log space
              box= c(FALSE,FALSE)
            )
          )
          
          sigmas=c(2,2,5,1,1,30,0.3,Inf)
          
          sigmas=list(
            sersic=list(
              xcen= numeric(2)+sigmas[1],
              ycen= numeric(2)+sigmas[2],
              mag= numeric(2)+sigmas[3],
              re= numeric(2)+sigmas[4],
              nser= numeric(2)+sigmas[5],
              ang= numeric(2)+sigmas[6],
              axrat= numeric(2)+sigmas[7],
              box= numeric(2)+sigmas[8]
            )
          )
          
          priors=profitMakePriors(modellist, sigmas, tolog, allowflat=TRUE)
          
          # The hard intervals should also be specified in log space if relevant:
          
          intervals=list(
            sersic=list(
              xcen=list(lim=c(main_src$xcen-20,main_src$xcen+20),lim=c(main_src$xcen-20,main_src$xcen+20)),
              ycen=list(lim=c(main_src$ycen-20,main_src$ycen+20),lim=c(main_src$ycen-20,main_src$ycen+20)),
              mag=list(lim=c(10,50),lim=c(15,50)),
              re=list(lim=c(1,100),lim=c(1,100)),
              nser=list(lim=c(0.5,10),lim=c(0.5,2)),
              ang=list(lim=c(-180,360),lim=c(-180,360)),
              axrat=list(lim=c(0.1,1),lim=c(0.1,1)),
              box=list(lim=c(-1,1),lim=c(-1,1))
            )
          )
          
          # constraints -------------------------------------
          #         constraints=function(modellist){
          #           if(modellist$sersic$re[1]>modellist$sersic$re[2]){
          #             modellist$sersic$re[1]=modellist$sersic$re[2]
          #           }
          #           return=modellist
          #         }
          
          #tempCL=profitOpenCLEnv()
          Data <- profitSetupData(image=image, mask=mask, sigma=sigma, segim=segim,
                                  psf=psf, modellist=modellist, tofit=tofit, tolog=tolog,
                                  priors=priors, intervals=intervals, magzero=ZP,
                                  algo.func='optim', verbose=FALSE, like.func = 't')
          
          ##########################
          ###### Initial Model #####
          ##########################
          
          png(file = paste(output_dir,'/init_model_BT',BTratio_BRe[[ii]][2],
                           '_BRe',BTratio_BRe[[ii]][1],'.png',sep = ""),width=8,height=2.5,units="in",res=300)
          profitLikeModel(parm=Data$init,Data=Data,makeplots=TRUE,whichcomponents=list(sersic="all"))
          dev.off()
          
          # } else if(nComp ==1) {
          #   ang <- c(main_src$ang)
          #   xcen <- main_src$xcen
          #   ycen <- main_src$ycen
          #   # mag <- -2.5*log10(c(main_src$flux*0.1,main_src$flux*0.9))-48.6
          #   # mag <- main_src$mag
          #   flux <- main_src$flux
          #   mag <- -2.5*log10(main_src$flux)+ZP
          #   re=sqrt(main_src$N50/(pi*main_src$axrat))
          #   axrat<- main_src$axrat
          #   
          #   modellist=list(
          #     sersic=list(
          #       xcen= xcen,
          #       ycen= ycen,
          #       mag= mag ,
          #       re= re,
          #       nser= 1.0000,
          #       ang= ang, 
          #       axrat= axrat, #min/maj: 1= o, 0= |
          #       box=0     
          #     )
          #   )
          #   
          #   tofit=list(
          #     sersic=list(
          #       xcen= TRUE, # fit for xcen of both disk and bulge.
          #       ycen= TRUE, # fit for ycen of both disk and bulge.
          #       mag= TRUE, # Fit for both
          #       re= TRUE, # Fit for both
          #       nser= TRUE, # Fit for bulge
          #       ang= TRUE, # Fit for disk
          #       axrat= TRUE, # Fit for disk
          #       box= FALSE # Fit for neither
          #     )
          #   )
          #   
          #   # What parameters should be fitted in log space:
          #   
          #   tolog=list(
          #     sersic=list(
          #       xcen= FALSE,
          #       ycen= FALSE,
          #       mag= FALSE,
          #       re= TRUE,    #re is best fit in log space
          #       nser= TRUE,  #nser is best fit in log space
          #       ang= FALSE,
          #       axrat= TRUE, #axrat is best fit in log space
          #       box= FALSE
          #     )
          #   )
          #   
          #   sigmas=c(2,2,5,1,1,30,0.3,Inf)
          #   
          #   sigmas=list(
          #     sersic=list(
          #       xcen= sigmas[1],
          #       ycen= sigmas[2],
          #       mag= sigmas[3],
          #       re= sigmas[4],
          #       nser= sigmas[5],
          #       ang= sigmas[6],
          #       axrat= sigmas[7],
          #       box= sigmas[8]
          #     )
          #   )
          #   
          #   priors=profitMakePriors(modellist, sigmas, tolog, allowflat=TRUE)
          #   
          #   # The hard intervals should also be specified in log space if relevant:
          #   
          #   intervals=list(
          #     sersic=list(
          #       xcen=list(lim=c(main_src$xcen-20,main_src$xcen+20)),
          #       ycen=list(lim=c(main_src$ycen-20,main_src$ycen+20)),
          #       mag=list(lim=c(10,30)),
          #       re=list(lim=c(1,100)),
          #       nser=list(lim=c(0.5,10)),
          #       ang=list(lim=c(-180,360)),
          #       axrat=list(lim=c(0.1,1)),
          #       box=list(lim=c(-1,1))
          #     )
          #   )
          #   
          #   #tempCL=profitOpenCLEnv()
          #   Data <- profitSetupData(image=image, mask=mask, sigma=sigma, segim=segim,
          #                           psf=psf, modellist=modellist, tofit=tofit, tolog=tolog,
          #                           priors=priors, intervals=intervals, magzero=ZP,
          #                           algo.func='optim', verbose=FALSE, like.func = 't')
          #   
          #   ##########################
          #   ###### Initial Model #####
          #   ##########################
          #   
          #   png(file = paste(output_dir,'initial_model.png',sep = "/"),width=8,height=2.5,units="in",res=300)
          #   profitLikeModel(parm=Data$init,Data=Data,makeplots=TRUE,whichcomponents=list(sersic=1))
          #   dev.off()
          #   
          # } # end of if on fit components.
          # 
          
          ##########################
          ###### Initial Model ellipse plot #####
          ##########################
          
          png(file = paste(output_dir,'/init_ellipse_BT',BTratio_BRe[[ii]][2],
                           '_BRe',BTratio_BRe[[ii]][1],'.png',sep = ""),width=7,height=5,units="in",res=300)
          if (nComp == 1){
            try(profitEllipsePlot(Data=Data,modellist=add_pseudo_bulge(modellist),pixscale=PIXELSCALE,FWHM=0.124,SBlim=26))
          } else if (nComp == 2){
            try(profitEllipsePlot(Data=Data,modellist=modellist,pixscale=PIXELSCALE,FWHM=0.124,SBlim=26))
          }
          dev.off()
          
          log7=" Initial Model:: DONE :D "

          # Initial values to be written in the output catalog #
          Init_mag1=Data$modellist$sersic$mag[1]
          if (nComp==2) Init_mag2=Data$modellist$sersic$mag[2] else Init_mag2=0 
          Init_re1=Data$modellist$sersic$re[1]
          if (nComp==2) Init_re2=Data$modellist$sersic$re[2]   else Init_re2=0
          Init_n1=Data$modellist$sersic$nser[1]
          if (nComp==2) Init_n2=Data$modellist$sersic$nser[2]  else Init_n2=0
          Init_ang1=Data$modellist$sersic$ang[1]
          if (nComp==2) Init_ang2=Data$modellist$sersic$ang[2] else Init_ang2=0
          Init_axrat1=Data$modellist$sersic$axrat[1]
          if (nComp==2) Init_axrat2=Data$modellist$sersic$axrat[2] else Init_axrat2=0          
          #######
          
          #########################
          ####### optimfit ########
          #########################
          #       if(verbose)cat("Running optimfit ......", '\n')
          #       pdf(file = paste(output_dir,'/optim_fit.pdf',sep = ""),height = 5)
          #       
          #       optimfit=optim(Data$init, profitLikeModel, method='BFGS',
          #                      Data=Data, control=list(fnscale=-1)) #,parscale=sigmas[which(unlist(tofit))]   method='L-BFGS-B'
          #       profitLikeModel(optimfit$par,Data,makeplots=TRUE,whichcomponents=list(sersic='all'))
          #       profitLikeModel(optimfit$par,Data,makeplots=TRUE,whichcomponents=list(sersic='all'),plotchisq = T)
          #       #lower=lowers[which(unlist(tofit))], upper=uppers[which(unlist(tofit))],
          #       modeloptim=profitRemakeModellist(parm=optimfit$par, Data=Data)$modellist
          #       if (nComps == 1){
          #         try(profitEllipsePlot(Data=Data,modellist=add_pseudo_bulge(modeloptim),pixscale=PIXELSCALE,FWHM=0.124,SBlim=26))
          #       } else if (nComps == 2){
          #         try(profitEllipsePlot(Data=Data,modellist=modeloptim,pixscale=PIXELSCALE,FWHM=0.124,SBlim=26))
          #       }
          #       dev.off()
          #       if(verbose)cat("optimfit finished ......", '\n')
          #       log8=" optimfit:: DONE :D "
          
          
          ##########  LA  #############
          if(verbose)cat("Running Laplace Approximation ......", '\n')
          #       pdf(file = paste(output_dir,'/LAfit.pdf',sep = ""), width = 10, height = 8) 
          Data$algo.func="LA"
          
          LAfit=LaplaceApproximation(profitLikeModel, parm=Data$init, Data=Data, Iterations=1e3,
                                     Method='LM', CovEst='Identity', sir=FALSE)
          #       profitLikeModel(LAfit$Summary1[,1],Data,makeplots=TRUE,whichcomponents=list(sersic='all'))
          #       profitLikeModel(LAfit$Summary1[,1],Data,makeplots=TRUE,whichcomponents=list(sersic='all'),plotchisq = T)
          modelLA=profitRemakeModellist(LAfit$Summary1[,1],Data$modellist,Data$tofit,Data$tolog)$modellist
          #       if (nComps == 1){
          #         try(profitEllipsePlot(Data=Data,modellist=add_pseudo_bulge(modelLA),pixscale=PIXELSCALE,FWHM=0.124,SBlim=26))
          #       } else if (nComps == 2){
          #         try(profitEllipsePlot(Data=Data,modellist=modelLA,pixscale=PIXELSCALE,FWHM=0.124,SBlim=26))
          #       }
          #       dev.off()
          log9=" Laplace Approximation:: DONE :D "
          
          ### Check for divergence
          
          
          if(verbose)cat("Testing LAFit for divergence....",'\n')
          log09=" Testing LAfit for divergence... :D "
          LAConverge = TRUE # IF converged=TRUE THEN use these as initial conditions. ELSE run LaplacesApproximation()
          for (n in seq(nComp)){
            if (is.element(modelLA$sersic$mag[n],Data$intervals$sersic$mag[[n]])){LAConverge=FALSE}
            if (is.element(modelLA$sersic$re[n],Data$intervals$sersic$re[[n]])){LAConverge=FALSE}
          }
          
          # If solution has converged, set the LAFit solution as the initial conditions
          if (LAConverge==TRUE){
            if(verbose)cat("LAfit solutions will be used as LD initial conditions.",'\n')
            log009="LA converged: LAfit solutions will be used as LD initial conditions."
            LD_InitVal = LAfit$Summary1[,1]
          } else {
            LD_InitVal=Data$init
            if(verbose)cat(" LAFit solution did not converge.",'\n')
            log009="LAfit solution did not converge. :("
          }
          
          #########################
          ######### MCMC ##########
          #########################
          #     
          #     is.model(profitLikeModel, Initial.Values=Data$init, Data=Data)
          #      LDfit=LaplacesDemon.hpc(profitLikeModel, Initial.Values=Data$init, Data=Data,
          #                              Iterations=1e4, Algorithm='CHARM', Thinning=1, 
          #                              Specs=list(alpha.star=0.44), Chains=10, CPUs=10)
          
          ##########  LD  #############
          ##########  Componentwise Hit-And-Run Metropolis (CHARM)  #############
          if(verbose)cat("Running FULL MCMC: Laplaces Demon ......", '\n')
          #       pdf(file = paste(output_dir,'/LDfit.pdf',sep = ""), width = 10, height = 8)
          Data$algo.func="LD"    
          LDfit=LaplacesDemon(profitLikeModel, Initial.Values=LD_InitVal, Data=Data,
                              Iterations=1e4, Algorithm='CHARM', Thinning=1, Specs=list(alpha.star=0.44))  #c(rep(0,14))
          
          png(file = paste(output_dir,'/LD_Post_all_BT',BTratio_BRe[[ii]][2],
                                  '_BRe',BTratio_BRe[[ii]][1],'.png',sep = ""),width=20,height=20,units="in",res=300)
          BestLD=magtri(LDfit$Posterior1, samples=1000, samptype='end')
          dev.off()
          #       system(paste('rm ',output_dir,'/LD_post_all.png', sep=""))
          png(file = paste(output_dir,'/LD_BT',BTratio_BRe[[ii]][2],
                           '_BRe',BTratio_BRe[[ii]][1],'.png',sep = ""),width=8,height=2.5,units="in",res=300)
          profitLikeModel(BestLD[,1],Data,makeplots=TRUE,whichcomponents=list(sersic='all'))
          dev.off()
          png(file = paste(output_dir,'/LD_chisq_BT',BTratio_BRe[[ii]][2],
                           '_BRe',BTratio_BRe[[ii]][1],'.png',sep = ""),width=10,height=7,units="in",res=300)
          profitLikeModel(BestLD[,1],Data,makeplots=TRUE,whichcomponents=list(sersic='all'),plotchisq = T)
          dev.off()
          
          modelLD=profitRemakeModellist(BestLD[,1],Data$modellist,Data$tofit,Data$tolog)$modellist
          
          png(file = paste(output_dir,'/LD_ellipse_BT',BTratio_BRe[[ii]][2],
                           '_BRe',BTratio_BRe[[ii]][1],'.png',sep = ""),width=10,height=8,units="in",res=300)
          if (nComp == 1){
            try(profitEllipsePlot(Data=Data,modellist=add_pseudo_bulge(modelLA),pixscale=PIXELSCALE,FWHM=0.124,SBlim=26))
          } else if (nComp == 2){
            try(profitEllipsePlot(Data=Data,modellist=modelLA,pixscale=PIXELSCALE,FWHM=0.124,SBlim=26))
          }
          dev.off()
          
          png(file = paste(output_dir,'/LD_Post1_BT',BTratio_BRe[[ii]][2],
                           '_BRe',BTratio_BRe[[ii]][1],'.png',sep = ""),width=10,height=10,units="in",res=300)
          magtri(LDfit$Posterior1[,3:ncol(LDfit$Posterior1)], samples=1000, samptype='end')
          dev.off()
          #       png(file = paste(output_dir,'LD_post2.png',sep = "/"),width=10,height=10,units="in",res=300)
          #       magtri(LDfit$Posterior2[,5:ncol(LDfit$Posterior1)], samples=1000, samptype='end')
          #       dev.off()
          # Consort(LDfit)   # print the info of LD fit, the time of running (in minute) included.
          log11=" Laplaces Demon:: DONE :D "
          
          #       MCMC_xcen1=modelLD$sersic$xcen[1];MCMC_xcen2=modelLD$sersic$xcen[2]
          #       MCMC_ycen1=modelLD$sersic$ycen[1];MCMC_ycen2=modelLD$sersic$ycen[2]
          MCMC_mag1=modelLD$sersic$mag[1]
          if (nComp==2) MCMC_mag2=modelLD$sersic$mag[2] else MCMC_mag2=0 
          MCMC_re1=modelLD$sersic$re[1]
          if (nComp==2) MCMC_re2=modelLD$sersic$re[2]   else MCMC_re2=0
          MCMC_n1=modelLD$sersic$nser[1]
          if (nComp==2) MCMC_n2=modelLD$sersic$nser[2]  else MCMC_n2=0
          MCMC_ang1=modelLD$sersic$ang[1]
          if (nComp==2) MCMC_ang2=modelLD$sersic$ang[2] else MCMC_ang2=0
          MCMC_axrat1=modelLD$sersic$axrat[1]
          if (nComp==2) MCMC_axrat2=modelLD$sersic$axrat[2] else MCMC_axrat2=0
          #       MCMC_box1=modelLD$sersic$box[1]
          #       if (nComp==2) MCMC_box2=modelLD$sersic$box[2] else MCMC_box2=0
          
          #############################################
          ######## Make models for all images save and Surface Brightness Profiles #########
          #############################################
          noise=matrix(rnorm(nrow(image)*nrow(image), mean=main_src$sky_mean, sd=main_src$skyRMS_mean), 
                       nrow(image), nrow(image))
          
          Init_model=profitRemakeModellist(Data$init,Data$modellist,Data$tofit,Data$tolog)$modellist
          
          Initial_model <- profitMakeModel(psf=psf, modellist=Init_model,dim = dim(image),
                                           whichcomponents = list(sersic="all"), magzero=ZP)
          
          LDmodel <- profitMakeModel(psf=psf, modellist=modelLD,dim = dim(image),
                                     whichcomponents = list(sersic="all"), magzero=ZP)
          
          LDmodel_noise=LDmodel$z+noise
          
          if (nComp == 2 ) {
            bulge_model <- profitMakeModel(psf=psf, modellist=modelLD,dim = dim(image),
                                           whichcomponents = list(sersic=1), magzero=ZP)
            disk_model <- profitMakeModel(psf=psf, modellist=modelLD,dim = dim(image),
                                          whichcomponents = list(sersic=2), magzero=ZP)
            
            disk_model_noise=disk_model$z+noise
            bulge_model_noise=bulge_model$z+noise
          } 
          
          
          
          # Save Image, Initial model, Fit model and Fit model+Random Nois ##############################
          cmap = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
          medimg = median(abs(image[image > 0 & Data$region]))
          maximg = max(abs(image[Data$region]))
          stretchscale = 1/medimg
          
          png(file = paste(output_dir,'/all_images_BT',BTratio_BRe[[ii]][2],
                           '_BRe',BTratio_BRe[[ii]][1],'.png',sep = ""),width=10.4,height=3,units="in",res=300)
          par(mfrow=c(1,4))
          #       layout(matrix(c(4,3,2,1), 1, 4, byrow = TRUE), widths = 1, heights = c(2,2,2,2,2), respect = FALSE)
          #       par(mar = c(3, 0, 4,2))  #c(bottom, left, top, right) 
          magimage(image, stretchscale = stretchscale, lo = -maximg, 
                   hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
                   xlab = "x/pix", ylab = "y/pix")
          mtext("Image", side=3)
          #       par(mar = c(3, 0, 4, 0))
          magimage(Initial_model, stretchscale = stretchscale, lo = -maximg, 
                   hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
                   xlab = "x/pix", ylab = "y/pix")
          mtext("Initial Model", side=3)   
          magimage(LDmodel, stretchscale = stretchscale, lo = -maximg, 
                   hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
                   xlab = "x/pix", ylab = "y/pix")
          mtext("Final Model", side=3)
          magimage(LDmodel_noise, stretchscale = stretchscale, lo = -maximg, 
                   hi = maximg, type = "num", zlim = c(0, 1), col = cmap, 
                   xlab = "x/pix", ylab = "y/pix")
          mtext("Final Model+Noise", side=3)
          dev.off()
          
          # Surface Brightness Profiles #######################################################
          # png(file = paste(output_dir,'/profile_MCMC.png',sep = ""), width = 10, height = 8, units="in", res=300)
          # #       pdf(file = paste(output_dir,'profile_MCMC.pdf',sep = "/"), width = 10, height = 8)
          # HST_mu(image =  image, main_source =  main_src, segim = segim ,model =  model_noise, modelOPlot = TRUE)
          # HST_mu(image=image,main_source=main_src,segim=segim,model=LDmodel_noise,nComp=nComp,
          #        col = 'green',modelOPlot=TRUE, legend=TRUE,main = "Surface Brightness Profile")
          # if (nComp == 2 ) {
          #   par(new=TRUE)
          #   HST_mu(image=image,main_source=main_src,segim=segim,model=bulge_model$z,
          #          modelOPlot=TRUE,col='red',legend = FALSE)
          #   par(new=TRUE)
          #   HST_mu(image=image,main_source=main_src,segim=segim,model=disk_model$z,
          #          modelOPlot=TRUE,col='blue',legend = FALSE)
          # }
          # 
          # dev.off()
          # log12=" Profile Plot:: DONE :D "
          # 
          # 
          ############################################
          ######### Write the output catalog #########
          ############################################
          
          t2 = (proc.time()-t1)/60
          elapsed_time=t2[3]
          
          LDmodel_path= paste(output_dir,'/LD_BT',BTratio_BRe[[ii]][2],
                                    '_BRe',BTratio_BRe[[ii]][1],'.png',sep = "") # to write in the last column of the output catalog
          
          cat_list <- data.frame(CATAID_out,z_out,look_back_t_out,RA_out,DEC_out,LAMBDAR_maj_ax_out,Ymag_out,stellar_mass_out,
                                 ProFound_R50_out,ProFound_R100_out,ProFound_maj_ax_out,ProFound_mag_out,edge_frac_out,asym_out,nComp,
                                 BTratio_BRe[[ii]][2],(1-BTratio_BRe[[ii]][2]),
                                 Init_mag1,Init_mag2,Init_re1,
                                 Init_re2,Init_n1,Init_n2,Init_ang1,Init_ang2,
                                 Init_axrat1,Init_axrat2,
                                 MCMC_mag1,MCMC_mag2,MCMC_re1,
                                 MCMC_re2,MCMC_n1,MCMC_n2,MCMC_ang1,MCMC_ang2,
                                 MCMC_axrat1,MCMC_axrat2,elapsed_time,LDmodel_path)
          
          write.table(cat_list, file = out_catalog, append = TRUE,
                      col.names = FALSE, row.names = FALSE,sep = ",")
          
          log13= " Output Catalogue:: DONE :D "
          log14=paste(' Elapsed time:: ',t2[3])
          log15=" ----------------------------- "
          cat(log15,'\n',log01,'\n',log1,'\n',log001,'\n',log2, '\n',log3,'\n',log4,'\n',log04,'\n',log004,'\n',log5,'\n',log6,'\n',
              log7,'\n',log8,'\n',log9,'\n',log09,'\n',log009,'\n',log11,'\n',log12,'\n',log13,'\n',
              log14,'\n', file=paste(wrk_dir,'/',logfile, sep = '' ), append=TRUE)
          
        } # End of loop on BTratio and BulgeRe.
        
    }, error=function(e){
      t2 = (proc.time()-t1)/60
      elapsed_time=t2[3]
      print("error")
      
      LDmodel_path= paste(output_dir,'/LD_BT',BTratio_BRe[[ii]][2],
                          '_BRe',BTratio_BRe[[ii]][1],'.png',sep = "") # to write in the last column of the output catalog
      
      err_val <- data.frame(CATAID_out,z_out,look_back_t_out,RA_out,DEC_out,LAMBDAR_maj_ax_out,Ymag_out,stellar_mass_out,
                            ProFound_R50_out,ProFound_R100_out,ProFound_maj_ax_out,ProFound_mag_out,edge_frac_out,asym_out,nComp,
                            BTratio_BRe[[ii]][2],(1-BTratio_BRe[[ii]][2]),
                            Init_mag1,Init_mag2,Init_re1,
                            Init_re2,Init_n1,Init_n2,Init_ang1,Init_ang2,
                            Init_axrat1,Init_axrat2,
                            MCMC_mag1,MCMC_mag2,MCMC_re1,
                            MCMC_re2,MCMC_n1,MCMC_n2,MCMC_ang1,MCMC_ang2,
                            MCMC_axrat1,MCMC_axrat2,elapsed_time,LDmodel_path)
      
      write.table(err_val, file = out_catalog, append = TRUE, 
                  col.names = FALSE, row.names = FALSE,sep = ",")
      
      log13= " Output Catalogue:: DONE :D "
      
      log0= '*** ERROR ***'
      log14=paste(' Elapsed time:: ',t2[3])
      log15=" ----------------------------- "
      
      cat(log15,'\n',log0,'\n',log01,'\n',log1,'\n',log001,'\n',log2, '\n',log3,'\n',log4,'\n',log04,'\n',log004,'\n',log5,'\n',log6,'\n',
          log7,'\n',log8,'\n',log9,'\n',log09,'\n',log009,'\n',log11,'\n',log12,'\n',log13,'\n',
          log14,'\n', file=paste(wrk_dir,'/',logfile, sep = '' ), append=TRUE)
      
    })    # End of tryCatch
  # }       # End of the loop on galaxies 
  
  header <- c("CATAID","z","look_back_t","RA","DEC","LAMBDAR_maj_ax[asec]","Ymag","stellar_mass[solar]",
              "ProFound_R50[asec]","ProFound_R100[asec]","ProFound_maj_ax[pixel]","ProFound_mag","ProFound_frac","ProFound_asym","nComp",
              "BTflux","DTflux",
              "Init_mag1","Init_mag2","Init_re1[pixel]",
              "Init_re2[pixel]","Init_n1","Init_n2","Init_ang1","Init_ang2",
              "Init_axrat1","Init_axrat2",
              "MCMC_mag1","MCMC_mag2","MCMC_re1[pixel]",
              "MCMC_re2[pixel]","MCMC_n1","MCMC_n2","MCMC_ang1","MCMC_ang2",
              "MCMC_axrat1","MCMC_axrat2","elapsed_time[min]","LDmodel_path")
  
  cat<-read.csv(file = out_catalog, header = FALSE)
  write.table(cat, file = out_catalog, append = FALSE,
              col.names = header, row.names = FALSE,sep = ",")
  
  ############################################
  ############ Save the workspace ############
  ############################################
  
  # workspaceFilename = paste(output_dir,"/WorkSpace.RData", sep='')
  # save.image(workspaceFilename)
  
  return("Job Finished")
  cat("**** Job Finished ****",'\n' , file=paste(wrk_dir,'/',logfile, sep = '' ), append=TRUE)
  
}         
# END


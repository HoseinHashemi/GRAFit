#' GRAFit: Master code for connecting and calling other modules.
#'
#' @description This code is considered as the core code of the GRAFit that connects all other parts and calls modules where needed.
#' @param wrk_dir Working directory, where the GRAFit outputs will be saved.
#' @param data_dir The directory where the imaging data (frames) is stored.
#' @param PSF_dir The directory where pre-generated PSFs are stored. If this is not provided the PSF can be put into the wrk_dir. Alternatively, GRAFit will generate PSF for each galaxy.
#' @param object_list A catalogue of objects you wish to fit, should include: CATAID, RA, DEC, etc.
#' @param threadMode specifies the thereading mode of the parallel mode of the GRAFit. 0: snow-like, suitable for laptops 1: using Rmpi, suitable for large clusters like supercomputers.
#' @param ncores Number of CPUs to be used in parallel mode. Default = 1
#' @param logfile Logical; Should a log file be generated. Default = \code{TRUE}
#' @param nComp Number of components to be used in fitting. Default = 2 for a bulge+disk model. Alternatively could be 1 for a single Sersic model.
#' @param optimMode The optimization mode. Default = 'MCMC'; using \code{LaplacesDemon} package. Alternatives for this argument is: 'LA' that uses \code{Laplace Approximation} and 'optim' that uses \code{optim} optimization.
#' @param LA_iteration The iterations for \code{Laplace Approximation} optimization if optimMode = 'LA'. Default = 1000
#' @param Optim_algo The algorithm to be used in \code{optim} optimization. Default = \code{"BFGS"}. Available algorithms are: \code{c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")}
#' @param MCMC_algo The MCMC algorithm to be used in \code{Laplaces Demon}. Default = \code{"CHARM"}. See the \code{\link[LaplacesDemon]{LaplacesDemon}} documentation for all available algorithms.
#' @param MCMC_iteration The iterations for \code{MCMC} optimization if optimMode = 'MCMC'. Default = 10000
#' @param ExpDisk Logical; Should the disk be fitted with a pure exponential profile (\code{Sersic index = 1}). Default = \code{FALSE}; for a free Sersic profile.
#' @param FreeBulge Logical; Should the bulge location (\code{x \& y}) also be fitted by GRAFit. Default = \code{FALSE}, i.e. the bulge position will be fixed to the disk position.
#' @param BulgeFreeness In unit pixels. How many pixels the bulge position is allowed to be free from the centre of disk. Default = 11.
#' @param Single_PSF If \code{PSF_dir = FALSE} should the PSF be generated as a single PSF on the position of the galaxy on the CCD. Default = \code{FALSE}, generate a stacked PSF on the location of galaxy on each of the raw ACS images.
#' @param PSF_diameter the diameter of the PSF in arcseconds. Default = 1"
#' @param like.func likelihhod function to be parsed to the ProFit: "t": t-distribution, "norm": normal distribution. See ProFit package for more details.
#' @param catalog_name The name of GRAFit's output catalogue.
#' @param keep_wrk_space Logical; Should the work space be saved. Default: \code{TRUE}.
#' @param verbose Logical; Verbose.
#' @param plot Logical; Should plots be generated or just the structural catalogue.
#' @param add_hdr Logical; Should a header be appended to the top of the catalogue. Default = \code{TRUE}. Usefull for some large runs on supercomputers large dataset is splitted into several chunks.
#' @param DoPriors Logical; Should a prior distribution be applied to the fitting process. Default = \code{TRUE}.
#' @param DoConstraits Logical; Should constraints be applied to the fitting process. Default = \code{TRUE}.
#' @return GRAFit's output is a structural catalogue (catalog_name.csv) as well as a folder named with the object's ID containing plots and figures of the initial and final models.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitDynamo_v2}}
#' @examples
#' GRAFitMaster(wrk_dir = '~/Desktop/wrk_dir/', data_dir = '~/Desktop/data_dir/', PSF_dir = '~/Desktop/PSF_dir/', threadMode = 0, ncores = 1, nComp= 2, optimMode = 'MCMC', object_list = 'object_list')
#' @export

# devtools::install_github("renkun-ken/formattable")
# install_github('asgr/magicaxis',force=T)
# install_github("ICRAR/ProFit")
# install_github('asgr/ProFound')
# source("http://bioconductor.org/biocLite.R")
# biocLite("EBImage")

GRAFitMaster <- function( wrk_dir = NULL, data_dir = NULL, GRAFitlib = NULL, PSF_dir = NULL, object_list = NULL,
                         threadMode = c(0,1) , ncores = 1, logfile = "logfile.txt",
                         nComp = 2, optimMode = "MCMC", LA_iteration = 1e3,
                         Optim_algo = 'BFGS', MCMC_algo = "CHARM",
                         MCMC_iteration = 1e4, ExpDisk = FALSE,
                         FreeBulge = FALSE, BulgeFreeness = 11, Single_PSF = FALSE, PSF_diameter = 1,
                         like.func = "t", catalog_name = 'MasterCat.csv', keep_wrk_space = FALSE,
                         verbose = TRUE, plot = FALSE, add_hdr = TRUE,
                         DoPriors = TRUE, DoConstraits = TRUE) {

  options(digits=10)

  # set working directory --------------
  if (is.null(wrk_dir)) setwd(getwd()) else setwd(wrk_dir)
  if (is.null(data_dir)) stop(" The data directory (dir in which science images are) should be provided. ")

  out_catalog = paste(wrk_dir,'/',catalog_name,sep = "")
  system(paste('mkdir ', wrk_dir,'/WrkSp', sep=''))
  WrkSp_dir = paste(wrk_dir,'/WrkSp', sep='')

  if (threadMode == 0) { # on pc, laptop or local clusters like MUNRO.
    registerDoParallel(cores = ncores)     # snow-like funtionality.
  } else if (threadMode == 1) { # using Rmpi on big clusters like Raijin or MAGNUS.
    registerDoSNOW(makeCluster(ncores))
  }

  # Check if all required libraries are installed. Install otherwise.
  if (!"devtools" %in% rownames(installed.packages())) install.packages('devtools', dependencies = T)
  if (!"knitr" %in% rownames(installed.packages())) install.packages('knitr', dependencies = T)
  if (!"ProFit" %in% rownames(installed.packages())) install_github("ICRAR/ProFit")
  if (!"FITSio" %in% rownames(installed.packages())) install.packages('FITSio', dependencies = T)
  if (!"ProFound" %in% rownames(installed.packages())) install_github('asgr/ProFound')
  if (!"astro" %in% rownames(installed.packages())) install.packages('astro', dependencies = T)
  if (!"celestial" %in% rownames(installed.packages())) install.packages('celestial', dependencies = T)
  if (!"magicaxis" %in% rownames(installed.packages())) install_github('asgr/magicaxis',force=T)
  if (!"foreach" %in% rownames(installed.packages())) install.packages('foreach', dependencies = T)
  if (!"doParallel" %in% rownames(installed.packages())) install.packages('doParallel', dependencies = T)
  if (!"LaplacesDemon" %in% rownames(installed.packages())) install.packages('LaplacesDemon', dependencies = T)
  if (!"RColorBrewer" %in% rownames(installed.packages())) install.packages('RColorBrewer', dependencies = T)
  if (!"Cairo" %in% rownames(installed.packages())) install.packages('Cairo', dependencies = T)
  # if (!"AllStarFit" %in% rownames(installed.packages())) install_github("taranu/AllStarFit")

  ######################################
  ########## PREPARE INPUT FOR PROFIT ##########
  ######################################
  error = "No error"
  foreach (j = 1:NROW(object_list)) %dopar% {     # Loop on all objects. ,.errorhandling="pass" , .combine = 'cbind', .multicombine = TRUE

    tryCatch({

      StartTime = timestamp()
      t1 <- proc.time()
      options(digits=10)
      assign("wrk_dir", wrk_dir, envir = .GlobalEnv)
      assign("data_dir", data_dir, envir = .GlobalEnv)
      assign("GRAFitlib", GRAFitlib, envir = .GlobalEnv)

      error = " ERROR: in taking object list column info."
      # if (missing(object_list)) stop("No object list. You should provide a valid object list including at least RA & DEC")
      if (is.null(object_list$D10CATAID)) {
        CATAID=rep(0,nrow(object_list)); cat("Warning: ID not found.")
        }  else CATAID=object_list$D10CATAID
      if (is.null(object_list$RA)) {
        RA=rep(0,nrow(object_list)); stop("No RA in the object_list")
        } else RA=object_list$RA
      if (is.null(object_list$DEC)) {
        DEC=rep(0,nrow(object_list)); stop("No DEC in the object_list")
        } else DEC=object_list$DEC
      if (is.null(object_list$ZBEST)) {
        z=rep(0,nrow(object_list)); cat("Warning: Redshift not found.")
        } else z=object_list$ZBEST
      if (is.null(object_list$lookbk_T)) {
        LBT=rep(0,nrow(object_list)); cat("Warning: Look back time not found.")
        } else LBT=object_list$lookbk_T
      if (is.null(object_list$SEMIMAJ)) {
        maj_ax=rep(0,nrow(object_list)); cat("Warning: Major axis not found.")
        } else maj_ax=object_list$SEMIMAJ
      if (is.null(object_list$YMAG)) {
        Ymag=rep(0,nrow(object_list)); cat("Warning: YMAG not found.")
        } else Ymag=object_list$YMAG
      if (is.null(object_list$STELLARMASS)) {
        stellar_mass=rep(0,nrow(object_list)); cat("Warning: Stellar mass not found.")
        } else stellar_mass=object_list$STELLARMASS

      ###############################################################
      log01=' CATAID :: NA :( '; log1=' RA & DEC:: NA :( '
      log001= 'Redshift:: NA :( '; log2=' Found galaxy in frame:: NA :( '
      log3=' Image cutout:: NA :( '; log4=' Segmentation Map:: NA :( '
      log04= " Main source found:: NA "; log004=" segmentation dilated:: NA :( "
      log5=' Sigma map:: NA :( '; log6=' PSF:: NA :( '; log7=' Initial Model:: NA :( '
      log8=' optimfit:: NA :( ';log9=' Laplace Approximation:: NA :( '
      log09=' Testing LA conversion:: NA  '; log009=' LAfit converged/diverged? NA :(';
      log11=' Laplaces Demon:: NA :( '; log12=' SBProfile Plot:: NA :( '
      log13=' Output Catalogue:: NA :( '; log14= ' Elapsed time:: NA :( '

      CATAID_out = NA; G10ID_out = NA; z_out = NA; look_back_t_out = NA; RA_out = NA
      DEC_out = NA; SEMIMAJ_out = NA; Ymag_out = NA; stellar_mass_out = NA
      ProFoundRA = NA; ProFoundDEC = NA; ProFound_R50_out = NA; ProFound_R90_out = NA; ProFound_R100_out = NA
      ProFound_semimaj_out = NA; ProFound_semimin_out = NA
      ProFound_mag_out = NA; ProFound_magErr_out = NA; edge_frac_out = NA; asym_out = NA; elapsed_time = NA

      optim_BTflux = NA; optim_DTflux = NA
      optim_xcen1 = NA; optim_xcen2 = NA
      optim_ycen1 = NA; optim_ycen2 = NA
      optim_mag1 = NA; optim_mag2 = NA
      optim_re1 = NA; optim_re2 = NA
      optim_n1 = NA; optim_n2 = NA
      optim_ang1 = NA; optim_ang2 = NA
      optim_axrat1 = NA; optim_axrat2 = NA
      optim_box1 = NA; optim_box2 = NA
      logLike = NA; dof = NA; AIC = NA;
      LML = NA; DIC1 = NA; DIC2 = NA; fitClass = NA

      SDxcen1 = NA; SDxcen2 = NA; SDycen1 = NA; SDycen2 = NA; SDmag1 = NA; SDmag2 = NA
      SDre1 = NA; SDre2 = NA; SDnser1 = NA; SDnser2 = NA; SDang1 = NA; SDang2 = NA
      SDaxrat1 = NA; SDaxrat2 = NA
      MCSExcen1 = NA; MCSExcen2 = NA; MCSEycen1 = NA; MCSEycen2 = NA; MCSEmag1 = NA
      MCSEmag2 = NA; MCSEre1 = NA; MCSEre2 = NA; MCSEnser1 = NA; MCSEnser2 = NA; MCSEang1 = NA
      MCSEang2 = NA; MCSEaxrat1 = NA; MCSEaxrat2 = NA

      ###############################################################

      ra_deg = RA[j]
      dec_deg = DEC[j]

      log01 = paste(' ID :: ',CATAID[j], sep = '')
      log1 = paste(' RA & DEC :: ',ra_deg,'+',dec_deg, sep = '')
      log001 = paste(' Redshift:: ', z[j])
      ra_form = as.numeric(format(ra_deg, digits = 8))
      dec_form = as.numeric(format(dec_deg, digits = 8))

      if (is.null(CATAID[j])) {
        system(paste('mkdir ', wrk_dir,'/',ra_form,'+',dec_form, sep=''))
        assign("output_dir", paste(wrk_dir,'/',ra_form,'+',dec_form, sep=''), envir = .GlobalEnv)
      } else if (!is.null(CATAID[j])) {
        system(paste('mkdir ', wrk_dir,'/', CATAID[j], sep=''))
        assign("output_dir", paste(wrk_dir,'/', CATAID[j], sep=''), envir = .GlobalEnv)
      }

      CATAID_out = CATAID[j]; G10ID_out = object_list$G10CATAID[j]; z_out = z[j]; look_back_t_out = LBT[j]; RA_out = ra_deg; DEC_out =dec_deg
      SEMIMAJ_out = maj_ax[j]; Ymag_out = Ymag[j]; stellar_mass_out = stellar_mass[j]

      # Find frame name.
      error = "ERROR in: Frame finding"
      frame_im_name = GRAFitFrameFinder_v2(GRAFitlib = GRAFitlib,
                                           data_dir = data_dir,
                                           target_loc = c(ra_deg, dec_deg))

      if(verbose) cat(' Found galaxy in frame:: ', frame_im_name, '\n')
      log2 = paste(' Found galaxy in frame:: ', frame_im_name)

      error = "ERROR in: reading frame"
      frame_im = readFITS( frame_im_name )
      frame_hdr = frame_im$hdr
      frame_im = frame_im$imDat

      # Read required info from header.
      PHOTFLAM = as.numeric( frame_hdr[which(frame_hdr == "PHOTFLAM")+1] )
      PHOTPLAM = as.numeric( frame_hdr[which(frame_hdr == "PHOTPLAM")+1] )
      PHOTZPT = as.numeric( frame_hdr[which(frame_hdr == "PHOTZPT")+1] )
      EXPTIME = as.numeric( frame_hdr[which(frame_hdr == "EXPTIME")+1] )
      EXPSTART = as.numeric( frame_hdr[which(frame_hdr == "EXPSTART")+1] )
      PHOTBW = as.numeric( frame_hdr[which(frame_hdr == "PHOTBW")+1] )
      CCDGAIN = as.numeric( frame_hdr[which(frame_hdr == "CCDGAIN")+1] )
      FILTER = as.character( frame_hdr[which(frame_hdr == "FILTER2")+1] )
      # pix_scale = profoundGetPixScale( frame_hdr )
      pix_scale = 0.03
      # DATE_OBS = as.character(frame_hdr[which(frame_hdr=="DATE-OBS")+1])

      # Make the electron map
      # Because images are in units of count/second then
      # we only need to multiply the pixel values by the exposure time! CCDGAIN is 1.

      gain_ADU <- EXPTIME * CCDGAIN    # gain_ADU = EXPTIME * CCD_GAIN (=1)
      frame_elec_map <- gain_ADU * frame_im
      #elec_mp=paste(wrk_dir,'/','input_fits','/','frame_elec_map.fits',sep="")
      #write.fits(frame_elec_map, file=elec_mp)
      # ZERO POINT ####################
      instr_ABMAG_ZPT = -2.5 * log10(PHOTFLAM) - (  5 * log10(PHOTPLAM) ) - 2.408
      #ZP = instr_ABMAG_ZPT + 2.5*log10(multiplier)
      ZP = instr_ABMAG_ZPT + 2.5*log10(gain_ADU)
      # ZP = 34.95672347

      #----Dynamic cut-out----------------------------------------------------
      error = "ERROR in: Dynamic cut out"
      dyn_cut = GRAFitDynamo_v2( im = frame_elec_map, loc = c(ra_deg, dec_deg), R90 = object_list$R90[j],
                                 StarPinPoint = FALSE, GRAFitlib = GRAFitlib, gain = 1,
                                 magzero = ZP, pix_scale = pix_scale, header = frame_hdr,
                                 ImPlot = plot, segPlot = FALSE, output_dir = output_dir )

      main_src = dyn_cut$main_src
      image = dyn_cut$sky_red_image
      segim = dyn_cut$segim
      mask = segim == 0
     # mask = dyn_cut$mask
      sep = main_src$sep  # Radial separation between the flux weighted centre (xcen & ycen) and the flux centre (xmax & ymax) [in asec]
      con = main_src$con  # concentration of the object (R50/R90)
      ProFoundRA = main_src$RAcen
      ProFoundDEC = main_src$Deccen

      main_src = cbind(CATAID_out, main_src)

      ProCathdr = colnames(main_src)
      colnames(main_src)[1] = "CATAID"

      ProCat = paste(wrk_dir,'/ProCat.csv',sep = "")

      if (!file.exists(ProCat)) {
        write.table( main_src, file = ProCat, append = TRUE,
                     col.names = TRUE, row.names = FALSE, sep = "," )
      } else {
        write.table( main_src, file = ProCat, append = TRUE,
                     col.names = FALSE, row.names = FALSE, sep = "," )
      }

      gc(verbose = F)

      main_srcSky = mean(dyn_cut$sky); main_srcSkyRMS = mean(dyn_cut$skyRMS)

      log3=" Image cutout:: DONE :D "
      log004=" dilation:: DONE :D "
      log4= " Source finding & Segmentation Map:: DONE :D "

      if (NROW(main_src) != 0) {
        log04= " Main source finding:: DONE :D "
      } else {
        log04= " No main source is found. :( "
      }

      ProFound_R50_out = main_src$R50; ProFound_R90_out = main_src$R90; ProFound_R100_out = main_src$R100
      ProFound_semimaj_out = main_src$semimaj*pix_scale; ProFound_semimin_out = main_src$semimin*pix_scale ;
      ProFound_mag_out = main_src$mag; ProFound_magErr_out = main_src$mag_err
      edge_frac_out = main_src$edge_frac; asym_out = main_src$asymm

      # myPlotSegIm(image, segim, sky = dyn_cut$seg$sky, main = main_src )

      # save the segmentation map as a png file.
#       jpeg(filename = paste(output_dir,'/seg.jpeg',sep = '')) #,width=7,height=7.6,units="in",res=150)
#       c <- profitSegimPlot(image,seg$segim)
#       dev.off()
      # if(verbose)cat('Seg Map saved:: DONE :D ', '\n', file=paste(wrk_dir,'/',logfile, sep = '' ), append=TRUE)

      # SIGMA map by Profound
      error = "ERROR in: sigma map"
      sigma = profoundMakeSigma(image, objects = dyn_cut$objects,
                            skyRMS = dyn_cut$skyRMS, image_units = 'elec', sky_units = 'elec',
                            read_units = 'elec', dark_units = 'elec', output_units = 'elec', plot = plot)

      log5= " Sigma map:: DONE :D "

      # Write sigma.fits
      # sig=paste(wrk_dir,'/','input_fits','/','sigma.fits',sep="")
      # write.fits(sigma, file=sig)
#
#       psf = GRAFitPSFgenerator(wrk_dir = wrk_dir, GRAFitlib = GRAFitlib, output_dir = output_dir, HST_Focus_val_file = data_dir,
#                                target_loc = c(ra_deg, dec_deg), header = frame_hdr, exBmV = object_list$EBV[j], verbose = verbose )

      error = "ERROR in: PSF generation."

      if ( "EBV" %in% colnames(object_list) ) {
        if ( !is.na(object_list$EBV[j]) ) {
          EBV = object_list$EBV[j]
        } else EBV = 0
      } else EBV = 0

    if ( !file.exists(paste(wrk_dir,"/","PSF.fits",sep = "")) & !file.exists(paste(PSF_dir, object_list$D10CATAID[j],"_PSF.fits",sep = "")) ) {
      if (Single_PSF) {
        error = "ERROR in: old PSF generation."
        PSF_name = "psf_1"
        psf1 = GRAFitPSFgenerator(wrk_dir = wrk_dir, GRAFitlib = GRAFitlib, output_dir = output_dir, PSF_name = PSF_name,
                                  HST_Focus_val_file = data_dir, header = frame_hdr, target_loc = c(object_list$X_1[j], object_list$Y_1[j]), loc_Unit = "xy",
                                  CCDCHIP = object_list$CCDCHIP_1[j], filter = "f814w", jitter = 3, ebmv = EBV,
                                  EXPSTART = EXPSTART, PSF_diameter= PSF_diameter, verbose = verbose )
        finalPSF = readFITS(paste(output_dir,"/",PSF_name,'_image00.fits', sep = ""))$imDat

      } else {
        error = "ERROR in: new PSF generation."
        tiny_SUB = 5
        PSF_name = "psf_1"
        psf1 = GRAFitPSFgenerator(wrk_dir = wrk_dir, GRAFitlib = GRAFitlib, output_dir = output_dir, PSF_name = PSF_name,
                                HST_Focus_val_file = data_dir, header = frame_hdr, target_loc = round(c(object_list$X_1[j], object_list$Y_1[j])), loc_Unit = "xy",
                                CCDCHIP = object_list$CCDCHIP_1[j], filter = "f814w", jitter = 3, ebmv = EBV,
                                EXPSTART = EXPSTART, PSF_diameter= PSF_diameter, SUB = tiny_SUB, verbose = verbose )
        psf1 = readFITS(paste(output_dir,"/",PSF_name,'_image00.fits', sep = ""))

        PSF_name = "psf_2"
        psf2 = GRAFitPSFgenerator(wrk_dir = wrk_dir, GRAFitlib = GRAFitlib, output_dir = output_dir, PSF_name = PSF_name,
                                  HST_Focus_val_file = data_dir, header = frame_hdr, target_loc = round(c(object_list$X_2[j], object_list$Y_2[j])), loc_Unit = "xy",
                                  CCDCHIP = object_list$CCDCHIP_2[j], filter = "f814w", jitter = 3, ebmv = EBV,
                                  EXPSTART = EXPSTART,PSF_diameter= PSF_diameter, SUB = tiny_SUB, verbose = verbose )
        psf2 = readFITS(paste(output_dir,"/",PSF_name,'_image00.fits', sep = ""))

        PSF_name = "psf_3"
        psf3 = GRAFitPSFgenerator(wrk_dir = wrk_dir, GRAFitlib = GRAFitlib, output_dir = output_dir, PSF_name = PSF_name,
                                  HST_Focus_val_file = data_dir, header = frame_hdr, target_loc = round(c(object_list$X_3[j], object_list$Y_3[j])), loc_Unit = "xy",
                                  CCDCHIP = object_list$CCDCHIP_3[j], filter = "f814w", jitter = 3, ebmv = EBV,
                                  EXPSTART = EXPSTART, PSF_diameter= PSF_diameter, SUB = tiny_SUB, verbose = verbose )
        psf3 = readFITS(paste(output_dir,"/",PSF_name,'_image00.fits', sep = ""))

        PSF_name = "psf_4"
        psf4 = GRAFitPSFgenerator(wrk_dir = wrk_dir, GRAFitlib = GRAFitlib, output_dir = output_dir, PSF_name = PSF_name,
                                  HST_Focus_val_file = data_dir, header = frame_hdr, target_loc = round(c(object_list$X_4[j], object_list$Y_4[j])), loc_Unit = "xy",
                                  CCDCHIP = object_list$CCDCHIP_4[j], filter = "f814w", jitter = 3, ebmv = EBV,
                                  EXPSTART = EXPSTART, PSF_diameter= PSF_diameter, SUB = tiny_SUB, verbose = verbose )
        psf4 = readFITS(paste(output_dir,"/",PSF_name,'_image00.fits', sep = ""))

        finalPSF <- GRAFitPSFMakeStack(psf1, psf2, psf3, psf4)
        error = finalPSF$error
        finalPSF <- finalPSF$finalPSF

        write.fits(finalPSF, file = paste(output_dir,"/","PSF.fits",sep = ""))

      }

      if(verbose) cat(" PSF generated :D" ,'\n')
      log6=" PSF :: DONE :D "
    } else if ( file.exists(paste(output_dir,"/","PSF.fits",sep = "")) ) {
      finalPSF <- readFITS(file = paste(output_dir,"/","PSF.fits",sep = ""))$imDat
      cat("PSF provided. It will be used! :D", '\n'); log6=" PSF :: DONE :D "
    } else if ( file.exists(paste(PSF_dir, object_list$D10CATAID[j],"_PSF.fits",sep = "")) ) {
      finalPSF <- readFITS(file = paste(PSF_dir, object_list$D10CATAID[j],"_PSF.fits", sep = ""))$imDat
      cat("PSF provided. It will be used! :D", '\n'); log6=" PSF :: DONE :D "
    } else if ( file.exists(paste(wrk_dir, "PSF.fits",sep = "")) ) {
      finalPSF <- readFITS(file = paste(wrk_dir, "PSF.fits",sep = ""))$imDat
      cat("PSF provided. It will be used! :D", '\n'); log6=" PSF :: DONE :D "
    }

print("Doing model setup")

      error = "ERROR in: model set up"
      modelSetUp = GRAFitInitModelSetupConst( image = image, GRAFitlib = GRAFitlib, mask = mask, sigma = sigma, segim = segim,
                             psf = finalPSF, output_dir = output_dir, main_src = main_src, pix_scale = pix_scale,
                             nComp = nComp, ZP = ZP, ExpDisk = ExpDisk, SBlim = 26,
                             FreeBulge = FreeBulge, BulgeFreeness = BulgeFreeness, like.func = like.func,
                             DoPriors = DoPriors, DoConstraits = DoConstraits )

      Data = modelSetUp$Data

      log7=" Initial Model:: DONE :D "

      ####### Optim Fit ########

      if (optimMode == "optim") {

        error = "ERROR in: model optim"
        optim = GRAFitOptim( output_dir = output_dir, GRAFitlib = GRAFitlib,
                             Model = profitLikeModel, Data = Data, verbose = verbose,
                             nComp = nComp, pixscale = pix_scale, FWHM = 0.09,
                             SBlim = 26, zeropoint = ZP, Algorithm = Optim_algo )

        optimfit <- optim$optimfit
        modeloptim <- optim$modeloptim

        if (nComp == 1) {
          fitClass = 7  # single-sersic
        } else {
          fitClass <- optim$SBprof$fitClass
        }


        log8 = " optimfit:: DONE :D "

        optim_xcen1 =  modeloptim$sersic$xcen[1]
        if (nComp == 2) optim_xcen2 = modeloptim$sersic$xcen[2] else optim_xcen2 = NA
          optim_ycen1 = modeloptim$sersic$ycen[1]
        if (nComp == 2) optim_ycen2 = modeloptim$sersic$ycen[2] else optim_ycen2 = NA
          optim_mag1 = modeloptim$sersic$mag[1]
        if (nComp == 2) optim_mag2 = modeloptim$sersic$mag[2] else optim_mag2 = NA
          optim_re1 = modeloptim$sersic$re[1]
        if (nComp == 2) optim_re2 = modeloptim$sersic$re[2]   else optim_re2 = NA
          optim_n1=modeloptim$sersic$nser[1]
        if (nComp == 2) optim_n2 = modeloptim$sersic$nser[2]  else optim_n2 = NA
          optim_ang1 = modeloptim$sersic$ang[1]
        if (nComp == 2) optim_ang2 = modeloptim$sersic$ang[2] else optim_ang2 = NA
          optim_axrat1 = modeloptim$sersic$axrat[1]
        if (nComp == 2) optim_axrat2 = modeloptim$sersic$axrat[2] else optim_axrat2 = NA

      } # end of if (optimMode == "optim")

      if (optimMode == "MCMC") {

        ##########  Laplace Approximation  #############
        # print("LA")
        # error = "ERROR in: LA"
        # LA = GRAFitLA(Model =  profitLikeModel,
        #               Initial.Values = Data$init,
        #               Data = Data,
        #               iteration = LA_iteration,
        #               verbose = verbose)

#         LAfit <- LA$LAfit
#         modelLA <- LA$modelLA
#
#         log9=" Laplace Approximation:: DONE :D "
#
#         ### Check for LA convergence.
#
#         # if(verbose) cat( "Testing LAFit for divergence....",'\n' )
#         # log09 = " Testing LAfit for divergence... :D "
#         # LAConverge = TRUE # IF LAConverged = TRUE then use that as initial values for LD ELSE use the very initial guesses from ProFound.
#         # for (n in seq(nComp)){
#         #   if (is.element(format(modelLA$sersic$mag[n], digits = 6), Data$intervals$sersic$mag[[n]])) { LAConverge = FALSE }
#         #   if (is.element(format(modelLA$sersic$re[n], digits = 6), Data$intervals$sersic$re[[n]])) { LAConverge = FALSE }
#         # }
#         #
#         # # If solution converged, set the LAFit solution as the initial values fro LD.
#         # if ( LAConverge == TRUE ){
#         #   if( verbose ) cat(" LA converged: LA solutions will be used as LD initial conditions. :)",'\n')
#         #   log009 = " LA converged: LA solutions will be used as LD initial conditions. :)"
#         #   LD_InitVal = LAfit$Summary1[,1]
#         # } else {
#         #   LD_InitVal = Data$init
#         #   if(verbose) cat(" LA solution did not converge. :(",'\n')
#         #   log009 = " LA solution did not converge. :("
#         # }
#
#         LD_InitVal = LAfit$Summary1[,1]
          LD_InitVal = NULL
        ######### MCMC ##########
        #########  Laplaces Demon  ##########
        print("LD")

        error = "ERROR in: LD"
        LD = GRAFitLD( output_dir = output_dir, GRAFitlib = GRAFitlib, Model = profitLikeModel, Initial.Values = LD_InitVal,
                  Data = Data, iteration = MCMC_iteration, verbose = verbose, nComp = nComp, zeropoint = ZP,
                  pixscale = pix_scale, FWHM = 0.09, SBlim = 26, Algorithm = MCMC_algo  )

        LDfit <- LD$LDfit
        modelLD <- LD$modelLD

        if (nComp == 1) {
          fitClass = 7  # single-sersic
        } else {
          fitClass <- LD$SBprof$fitClass
        }

        logLike = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "LL", 1]
        dof = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "dof", 1]

        AIC = -2*logLike+2*dof
        LML = LDfit$LML
        if (is.na(LML)) LML = LML(LL=LDfit$Deviance*(-1/2), method="HME")$LML   # Taken from the examples of the LML func. Double Check !!
        DIC1 = LDfit$DIC1[3]
        DIC2 = LDfit$DIC2[3]
        # ERRORs
        SD = matrix(LD$LDfit$Summary1[, 2], 1, length(LD$LDfit$Summary1[, 2]))   # standard deviations
        MCSE = matrix(LD$LDfit$Summary1[, 3], 1, length(LD$LDfit$Summary1[, 3])) # Monte Carlo Standard Error
        log11 = ' Laplaces Demon :: DONE :D '

        SDxcen1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.xcen1", 2]; if (length(SDxcen1)==0) SDxcen1 = NA
        SDxcen2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.xcen2", 2]; if (length(SDxcen2)==0) SDxcen2 = NA
        SDycen1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.ycen1", 2]; if (length(SDycen1)==0) SDycen1 = NA
        SDycen2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.ycen2", 2]; if (length(SDycen2)==0) SDycen2 = NA
        SDmag1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.mag1", 2]; if (length(SDmag1)==0) SDmag1 = NA
        SDmag2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.mag2", 2]; if (length(SDmag2)==0) SDmag2 = NA
        SDre1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.re1", 2]; if (length(SDre1)==0) SDre1 = NA
        SDre2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.re2", 2]; if (length(SDre2)==0) SDre2 = NA
        SDnser1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.nser1", 2]; if (length(SDnser1)==0) SDnser1 = NA
        SDnser2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.nser2", 2]; if (length(SDnser2)==0) SDnser2 = NA
        SDang1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.ang1", 2]; if (length(SDang1)==0) SDang1 = NA
        SDang2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.ang2", 2]; if (length(SDang2)==0) SDang2 = NA
        SDaxrat1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.axrat1", 2]; if (length(SDaxrat1)==0) SDaxrat1 = NA
        SDaxrat2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.axrat2", 2]; if (length(SDaxrat2)==0) SDaxrat2 = NA

        MCSExcen1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.xcen1", 3]; if (length(MCSExcen1)==0) MCSExcen1 = NA
        MCSExcen2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.xcen2", 3]; if (length(MCSExcen2)==0) MCSExcen2 = NA
        MCSEycen1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.ycen1", 3]; if (length(MCSEycen1)==0) MCSEycen1 = NA
        MCSEycen2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.ycen2", 3]; if (length(MCSEycen2)==0) MCSEycen2 = NA
        MCSEmag1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.mag1", 3]; if (length(MCSEmag1)==0) MCSEmag1 = NA
        MCSEmag2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.mag2", 3]; if (length(MCSEmag2)==0) MCSEmag2 = NA
        MCSEre1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.re1", 3]; if (length(MCSEre1)==0) MCSEre1 = NA
        MCSEre2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.re2", 3]; if (length(MCSEre2)==0) MCSEre2 = NA
        MCSEnser1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.nser1", 3]; if (length(MCSEnser1)==0) MCSEnser1 = NA
        MCSEnser2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.nser2", 3]; if (length(MCSEnser2)==0) MCSEnser2 = NA
        MCSEang1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.ang1", 3]; if (length(MCSEang1)==0) MCSEang1 = NA
        MCSEang2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.ang2", 3]; if (length(MCSEang2)==0) MCSEang2 = NA
        MCSEaxrat1 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.axrat1", 3]; if (length(MCSEaxrat1)==0) MCSEaxrat1 = NA
        MCSEaxrat2 = LD$LDfit$Summary1[row.names(LD$LDfit$Summary1) == "sersic.axrat2", 3]; if (length(MCSEaxrat2)==0) MCSEaxrat2 = NA

        optim_xcen1 = modelLD$sersic$xcen[1]
        if (nComp == 2 & FreeBulge) {optim_xcen2 = modelLD$sersic$xcen[2]} else {optim_xcen2 = NA}
        optim_ycen1 = modelLD$sersic$ycen[1]
        if (nComp == 2 & FreeBulge) {optim_ycen2 = modelLD$sersic$ycen[2]} else {optim_ycen2 = NA}
        optim_mag1 = modelLD$sersic$mag[1]
        if (nComp == 2) {optim_mag2 = modelLD$sersic$mag[2]} else {optim_mag2 = NA}
        optim_re1 = modelLD$sersic$re[1]
        if (nComp == 2) {optim_re2 = modelLD$sersic$re[2]}   else {optim_re2 = NA}
        optim_n1 = modelLD$sersic$nser[1]
        if (nComp == 2) {optim_n2 = modelLD$sersic$nser[2]}  else {optim_n2 = NA}
        optim_ang1 = modelLD$sersic$ang[1]
        if (nComp == 2) {optim_ang2 = modelLD$sersic$ang[2]} else {optim_ang2 = NA}
        optim_axrat1 = modelLD$sersic$axrat[1]
        if (nComp == 2) {optim_axrat2 = modelLD$sersic$axrat[2]} else {optim_axrat2 = NA}
        # optim_box1 = modelLD$sersic$box[1]
        # if (nComp == 2) optim_box2 = modelLD$sersic$box[2] else optim_box2 = NA

      } # end of if (optimMode == "MCMC")

      # Bulge/Total & Disk/Total flux ratio.
      error = "ERROR in: B/T flux"
      if (nComp == 2 ) {
        flux1 = profoundMag2Flux( optim_mag1, magzero = ZP )
        flux2 = profoundMag2Flux( optim_mag2, magzero = ZP )
        tot_flux = flux1 + flux2
        BTflux = flux1/tot_flux
        # DTflux = flux2/tot_flux
      } else {BTflux = NA}

      t01 = (proc.time()-t1)/3600
      elapsed_time = t01[[3]] # the time upto this point.
      ###############################################################
      error = "ERROR in: fit_stat plotting"
      # png(file = paste(output_dir,"/","fit_stat.png",sep = ""), width = 21, height = 12, units = "cm", res = 200 )
      png( 1700, 1700, file = paste(output_dir,"/","ModelStat.png",sep = ""), res = 200)
      # par( mfcol = c(1,2) )
        layout(matrix(c(2,3,1,1), 2, 2 , byrow = F), widths = c(2,2,2), heights = c(2.27,2.27,2), respect = T)
        ylim = 13; xlim = 40
        plot(0, type='n', xaxt = 'n', yaxt = 'n', ann = FALSE, ylim = c(0, ylim), xlim = c(0, xlim))
        text(0.6, ylim, paste("ID = ", CATAID_out), adj = 0, cex = 1.3)
        text(0.6, ylim-1, paste("RA = ", format(RA_out, digits = 7) ), adj = 0, cex = 1.3)
        text(0.6, ylim-2, paste("DEC = ", format(DEC_out, digits = 5)), adj = 0, cex = 1.3)
        text(0.6, ylim-3, paste(object_list$ZTYPE[j], " = ", format(z_out, digits = 5)), adj = 0, cex = 1.3)
        text(0.6, ylim-4, paste("Total Mag = ", format(main_src$mag, digits = 5)), adj = 0, cex = 1.3)
        text(0.6, ylim-5, paste("log10(Stellar Mass) = ", object_list$STELLARMASS[j]), adj = 0, cex = 1.3)
        text(0.6, ylim-6, paste("B/T flux = ", format(BTflux, digits = 2)), adj = 0, cex = 1.3)

        # text(.89,0.2, paste("D/T flux = ", format(DTflux, digits = 2)))
        text(0.6,  ylim-7, paste("mag: ", "B = ",format(optim_mag1, digits = 4), " D = ",
                               format(optim_mag2, digits = 4)), adj = 0, cex = 1.3)
        # text(.89,-0.21, paste("D mag = ", format(optim_mag2, digits = 4)))
        text(0.6,  ylim-8, paste("Re/asec: ", "B = ", pix_scale*as.numeric(format(optim_re1, digits = 4))), adj = 0, cex = 1.3)
        text(0.6,  ylim-8.5, paste("                D = ", pix_scale*as.numeric(format(optim_re2, digits = 4))), adj = 0, cex = 1.3)
        # text(.93, -0.61, paste("D re = ", pix_scale*as.numeric(format(optim_re2, digits = 4)), "asec"))
        text(0.6,  ylim-9, paste("n: ", "B = ", format(optim_n1, digits = 3), "     D = ",
                               format(optim_n2, digits = 3)), adj = 0, cex = 1.3)
        text(0.6, ylim-10, paste("Fit Class = ", fitClass), adj = 0, cex = 1.3)
        # text(.85,-1.01, paste("D n1 = ", format(optim_n2, digits = 3)))
        text(0.6, ylim-11, paste("LL = ", format(logLike, digits = 5)), adj = 0, cex = 1.3)
        text(0.6, ylim-12, paste("DIC1 = ", format(DIC1, digits = 6)), adj = 0, cex = 1.3)
        text(0.6, ylim-13, paste("t = ", format(elapsed_time, digits = 5), "h"), adj = 0, cex = 1.3)

        magimage(dyn_cut$orig_image,
                 xlab = "x/pix", ylab = "y/pix",
                 cex.lab = 1.5,
                 cex.axis = 1.5,
                 mgp = c(2.2,1.2,0),
                 mtline = 3.5)
        title("Original")
        magimage(dyn_cut$sky_red_image,
                 xlab = "x/pix", ylab = "y/pix",
                 cex.lab = 1.5,
                 cex.axis = 1.5,
                 mgp = c(2.2,1.2,0),
                 mtline = 3.5)
        title("Background Subtracted")

      dev.off()

      #############################################
      ######## Make models for all images output and Surface Brightness Profiles #########
      #############################################

      error = "ERROR in: model generation for surface brightness profile"
      noise = matrix( rnorm( nrow(image) * nrow(image), mean = main_srcSky, sd = main_srcSkyRMS),
                   nrow(image), nrow(image))

      Init_model = profitRemakeModellist(Data$init, Data$modellist, Data$tofit, Data$tolog)$modellist

      Initial_model <- profitMakeModel(psf = finalPSF, modellist = Init_model, dim = dim(image),
                                      whichcomponents = list(sersic = "all"), magzero = ZP)

      if (optimMode == "optim") {
        final_model <- profitMakeModel(psf = finalPSF, modellist = modeloptim, dim = dim(image),
                                   whichcomponents = list(sersic = "all"), magzero = ZP)
      } else if (optimMode == "MCMC") {
        final_model <- profitMakeModel(psf = finalPSF, modellist = modelLD, dim = dim(image),
                                whichcomponents = list(sersic = "all"), magzero = ZP)
      }

      model_noise = final_model$z + noise

      if (nComp == 2 ) {
        if (optimMode == "optim") {
          bulge_model <- profitMakeModel( psf = finalPSF, modellist = modeloptim, dim = dim(image),
                                         whichcomponents = list(sersic = 1), magzero = ZP )
          disk_model <- profitMakeModel( psf = finalPSF, modellist = modeloptim, dim = dim(image),
                                        whichcomponents = list(sersic = 2), magzero = ZP )
        } else if (optimMode == "MCMC") {
          bulge_model <- profitMakeModel( psf = finalPSF, modellist = modelLD, dim = dim(image),
                                         whichcomponents = list(sersic = 1), magzero = ZP )
          disk_model <- profitMakeModel( psf = finalPSF, modellist = modelLD, dim = dim(image),
                                        whichcomponents = list(sersic = 2), magzero = ZP )
        }

        disk_model_noise = disk_model$z + noise
        bulge_model_noise = bulge_model$z + noise
      }

      # Surface Brightness Profiles #######################################################
      error = "ERROR in: surface brightness profile plotting"
      par( mfcol = c(1,1) )
      # png(file = paste(output_dir,'/SBprofile.png', sep = ""), width = 10, height = 8, units="in", res=300)
      png( 1000, 800, file = paste(output_dir,'/SBprofile.png', sep = ""), res = 100 )

        GRAFitSBprofile( image = image, main_source = main_src, segim = segim, model = model_noise, zeropoint = ZP,
                        comp = "bd", centerPos = c( optim_xcen1, optim_ycen1 ),
                        col = 'green', modelPlot = TRUE, legend = TRUE, title = "Surface Brightness Profile")

      invisible(
        if (nComp == 2 ) {
          par( new = TRUE )
          GRAFitSBprofile( image = image, main_source = main_src, segim = segim,
                          model = disk_model$z, zeropoint = ZP,
                          comp = "d", centerPos = c( optim_xcen1, optim_ycen1 ),
                          modelPlot = TRUE, col= 'blue', legend = FALSE )
          par( new = TRUE )
          GRAFitSBprofile( image = image, main_source = main_src, segim = segim,
                          model = bulge_model$z, zeropoint = ZP,
                          comp = "b", centerPos = c( optim_xcen1, optim_ycen1 ),
                          modelPlot = TRUE, col = 'red', legend = FALSE )
        }
      ) # end of invisible
      dev.off()

      log12=' SBProfile Plot:: DONE :D '


      t2 = (proc.time()-t1)/3600
      elapsed_time = t2[[3]]

      ############################################
      ######### Write the output catalog #########
      ############################################

      if (optimMode == "optim") {
        model_path=paste(output_dir,'/optim.png',sep='')
      } else model_path=paste(output_dir,'/LD.png',sep='')

        error = "ERROR in: writing out the master catalogue"

        cat_list <- data.frame( CATAID_out, G10ID_out, z_out, look_back_t_out, RA_out, DEC_out,
                                SEMIMAJ_out, Ymag_out, stellar_mass_out,
                                ProFoundRA = ProFoundRA, ProFoundDEC = ProFoundRA,
                                ProFound_R50_out, ProFound_R90_out, ProFound_R100_out, ProFound_semimaj_out, ProFound_semimin_out,
                                ProFound_mag_out, ProFound_magErr_out, edge_frac_out, asym_out, sep, con,
                                nComp, BTflux,
                                optim_xcen1, SDxcen1, MCSExcen1, optim_xcen2, SDxcen2, MCSExcen2,
                                optim_ycen1, SDycen1, MCSEycen1, optim_ycen2, SDycen2, MCSEycen2,
                                optim_mag1, SDmag1, MCSEmag1, optim_mag2, SDmag2, MCSEmag2,
                                optim_re1, SDre1, MCSEre1, optim_re2, SDre2, MCSEre2,
                                optim_n1, SDnser1, MCSEnser1, optim_n2, SDnser2, MCSEnser2,
                                optim_ang1, SDang1, MCSEang1, optim_ang2, SDang2, MCSEang2,
                                optim_axrat1, SDaxrat1, MCSEaxrat1, optim_axrat2, SDaxrat2, MCSEaxrat2,
                                logLike, dof, AIC, DIC1, DIC2, LML, fitClass, elapsed_time, model_path )


        write.table( cat_list, file = out_catalog, append = TRUE,
                    col.names = FALSE, row.names = FALSE, sep = ",")

        log13= " Output Catalogue:: DONE :D "
        log14=paste(' Elapsed time:: ',t2[3])
        log15=" ----------------------------- "
        cat(log15,'\n', StartTime,'\n',log01,'\n',log1,'\n',log001,'\n',log2, '\n',
            log3,'\n',log4,'\n',log04,'\n',log004,'\n',log5,'\n',log6,'\n',
            log7,'\n',log8,'\n',log9,'\n',log09,'\n',log009,'\n',log11,'\n',
            log12,'\n',log13,'\n', log14,'\n', file=paste(wrk_dir,'/',logfile, sep = '' ), append=TRUE)

        ############################################
        ############ Save the workspace ############
        ############################################
        if (keep_wrk_space) {
          workspaceFilename = paste( WrkSp_dir,"/",CATAID_out, ".RData", sep='' )
          save.image(workspaceFilename)
        }

    }, error=function(e){
      t2 = (proc.time()-t1)/60
      elapsed_time=t2[3]
      print(error)

      if (optimMode == "optim") {
        model_path=paste(output_dir,'/optim.png',sep='')
      } else model_path=paste(output_dir,'/LD.png',sep='')

      cat_list <- data.frame( CATAID_out, G10ID_out, z_out, look_back_t_out, RA_out, DEC_out,
                             SEMIMAJ_out, Ymag_out, stellar_mass_out,
                             ProFoundRA = ProFoundRA, ProFoundDEC = ProFoundDEC,
                             ProFound_R50_out, ProFound_R90_out, ProFound_R100_out, ProFound_semimaj_out, ProFound_semimin_out,
                             ProFound_mag_out, ProFound_magErr_out, edge_frac_out, asym_out, sep, con,
                             nComp, BTflux,
                             optim_xcen1, SDxcen1, MCSExcen1, optim_xcen2, SDxcen2, MCSExcen2,
                             optim_ycen1, SDycen1, MCSEycen1, optim_ycen2, SDycen2, MCSEycen2,
                             optim_mag1, SDmag1, MCSEmag1, optim_mag2, SDmag2, MCSEmag2,
                             optim_re1, SDre1, MCSEre1, optim_re2, SDre2, MCSEre2,
                             optim_n1, SDnser1, MCSEnser1, optim_n2, SDnser2, MCSEnser2,
                             optim_ang1, SDang1, MCSEang1, optim_ang2, SDang2, MCSEang2,
                             optim_axrat1, SDaxrat1, MCSEaxrat1, optim_axrat2, SDaxrat2, MCSEaxrat2,
                             logLike, dof, AIC, DIC1, DIC2, LML, fitClass, elapsed_time, model_path )

      write.table(cat_list, file = out_catalog, append = TRUE,
                    col.names = FALSE, row.names = FALSE, sep = ",")

      log13 = " Output Catalogue:: DONE :D "

        log0 = '*** ERROR ***'
        log14 = paste(' Elapsed time:: ', t2[3])
        log15 = " ----------------------------- "

        cat(log15,'\n', StartTime,'\n',log01,'\n',log1,'\n',log001,'\n',log2, '\n',
            log3,'\n',log4,'\n',log04,'\n',log004,'\n',log5,'\n',log6,'\n',
            log7,'\n',log8,'\n',log9,'\n',log09,'\n',log009,'\n',log11,'\n',
            log12,'\n',log13,'\n', log14,'\n', file=paste(wrk_dir,'/',logfile, sep = '' ), append=TRUE)

        ############################################
        ############ Save the workspace ############
        ############################################
        if (keep_wrk_space) {
          workspaceFilename = paste( WrkSp_dir,"/",CATAID_out, ".RData", sep='' )
          save.image( workspaceFilename )
        }


    })    # End of tryCatch
  }       # End of the loop on galaxies

  if (add_hdr) {

    GRAFitAddhdr(out_catalog = out_catalog)

     } # end of if(add_hdr)

  # splitNo = round(NROW(object_list)/50)
  # displayCom <- paste("/Users/22111305/anaconda/bin/python ", GRAFitlib, "/DisplayOutputs.py -D ", wrk_dir,
  #                     " --outputfile ", wrk_dir, "/Display.html",
  #                     " --run GRAFit", " --components ", nComp," --catFile ", wrk_dir,
  #                     "/catalog.csv", " --number ", splitNo  , sep = '' )
  # system(displayCom)

  cat("*** Job Finished ***", '\n')
  # cat("**** Job Finished ****",'\n' , file=paste(wrk_dir,'/',logfile, sep = '' ), append=TRUE)

  if (threadMode == 1) {
    stopCluster(cl)
    # mpi.exit()
  }

}
# END




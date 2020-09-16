GRAFitInitModelSetup <- function( output_dir = NULL, GRAFitlib = GRAFitlib, main_src = NULL, nComp = 2, 
                                  ExpDisk = FALSE, freeBulge = FALSE, BulgeFreeness = 3, BT_flux_ratio = 0.2, ZP = NULL, 
                                  image = image, mask = mask, sigma = sigma, pix_scale = pix_scale, 
                                  segim = segim, psf = psf, like.func = "t" )  {
  
  source(paste(GRAFitlib,'/GRAFitAddFakeBulge.R', sep=''))
  source(paste(GRAFitlib,'/GRAFitEllipsePlot.R',sep=''))
  
  if (missing(output_dir)) output_dir = getwd()
  
  if (nComp == 2) {
    ang <- c(0,main_src$ang)
    xcen <- main_src$xcen
    ycen <- main_src$ycen
    # xcen <- main_src$xmax
    # ycen <- main_src$ymax
    # mag <- -2.5*log10(c(main_src$flux*0.1,main_src$flux*0.9))-48.6
    # mag <- main_src$mag
    # flux <- main_src$flux
    mag <- -2.5 * log10( c( main_src$flux * BT_flux_ratio, main_src$flux * (1 - BT_flux_ratio) ) ) + ZP
    re = main_src$R50/pix_scale   # re = sqrt( main_src$N50 / (pi * main_src$axrat) ) [pixel]
    axrat <- c( 1, ( main_src$axrat ) )
    
    modellist = list(
      sersic = list(
        xcen = rep(xcen,2),
        ycen = rep(ycen,2),
        mag = mag ,
        re = c(0.2*re,re),
        nser = c(4., 1.),
        ang = ang, 
        axrat = axrat, #min/maj: 1= o, 0= |
        box = c(0., 0.)     
      )
    )
    
    if (ExpDisk) {
      if (freeBulge) {
        tofit = list(
          sersic = list(
            xcen = c(TRUE,TRUE),   # Fit for xcen of bulge and use it for disk as well.
            ycen = c(TRUE,TRUE),   # Fit for ycen of bulge and use it for disk as well.
            mag = c(TRUE,TRUE),  # Fit for both
            re = c(TRUE,TRUE),   # Fit for both
            nser = c(TRUE,FALSE), # Fit for bulge only
            ang = c(FALSE,TRUE), # Fit for disk
            axrat = c(FALSE,TRUE), # Fit for disk
            box = c(FALSE,FALSE) # Fit for neither
          )
        )  
      } else {
        tofit = list(
          sersic = list(
            xcen = c(TRUE,NA),   # Fit for xcen of bulge and use it for disk as well.
            ycen = c(TRUE,NA),   # Fit for ycen of bulge and use it for disk as well.
            mag = c(TRUE,TRUE),  # Fit for both
            re = c(TRUE,TRUE),   # Fit for both
            nser = c(TRUE,FALSE), # Fit for bulge only
            ang = c(FALSE,TRUE), # Fit for disk
            axrat = c(FALSE,TRUE), # Fit for disk
            box = c(FALSE,FALSE) # Fit for neither
          )
        )
      }
      
    } else {
      if (freeBulge) {
        tofit = list(
          sersic = list(
            xcen = c(TRUE,TRUE),   # Keep the centre of disk and bulge on top of each other.
            ycen = c(TRUE,TRUE),   # Keep the centre of disk and bulge on top of each other.
            mag = c(TRUE,TRUE),    # Fit for both
            re = c(TRUE,TRUE),     # Fit for both
            nser = c(TRUE,TRUE),   # Fit for both
            ang = c(FALSE,TRUE),   # Fit for disk
            axrat = c(FALSE,TRUE), # Fit for disk
            box = c(FALSE,FALSE)   # Fit for neither
          )
        )    
      } else {
        tofit = list(
          sersic = list(
            xcen = c(TRUE,NA),   # Keep the centre of disk and bulge on top of each other.
            ycen = c(TRUE,NA),   # Keep the centre of disk and bulge on top of each other.
            mag = c(TRUE,TRUE),    # Fit for both
            re = c(TRUE,TRUE),     # Fit for both
            nser = c(TRUE,TRUE),   # Fit for both
            ang = c(FALSE,TRUE),   # Fit for disk
            axrat = c(FALSE,TRUE), # Fit for disk
            box = c(FALSE,FALSE)   # Fit for neither
          )
        )    
      }
      
    }
    
    # What parameters should be fitted in log space:
    
    tolog = list(
      sersic = list(
        xcen = c(FALSE,FALSE),
        ycen = c(FALSE,FALSE),
        mag = c(FALSE,FALSE),
        re = c(TRUE,TRUE),    #re is best fit in log space
        nser = c(TRUE,TRUE),  #nser is best fit in log space
        ang = c(FALSE,FALSE),
        axrat = c(TRUE,TRUE), #axrat is best fit in log space
        box = c(FALSE,FALSE)
      )
    )
    
    sigmas = c(2,2,5,1,1,30,0.3,Inf)
    
    sigmas = list(
      sersic = list(
        xcen = numeric(2)+sigmas[1],
        ycen = numeric(2)+sigmas[2],
        mag = numeric(2)+sigmas[3],
        re = numeric(2)+sigmas[4],
        nser = numeric(2)+sigmas[5],
        ang = numeric(2)+sigmas[6],
        axrat = numeric(2)+sigmas[7],
        box = numeric(2)+sigmas[8]
      )
    )
    
    priors = profitMakePriors(modellist = modellist, sigmas = sigmas, 
                              tolog = tolog, tofit = tofit, allowflat = TRUE)
    
    # The hard intervals should also be specified in log space if relevant:
    freeness = BulgeFreeness -1
    freeness = freeness / 2
    
      intervals = list(
        sersic = list(
          xcen = list(lim = c(xcen - freeness, xcen + freeness), lim = c(xcen - 10.0, xcen + 10.0)),
          ycen = list(lim = c(ycen - freeness, ycen + freeness), lim = c(ycen - 10.0, ycen + 10.0)),
          mag = list(lim = c(10.0,40.0), lim = c(10.0,40.0)),
          re = list(lim = c(1.0,dim(image)[1]), lim = c(1.0,dim(image)[1])),
          nser = list(lim = c(0.5,20.0), lim = c(0.5,1.5)),
          ang = list(lim = c(-180.,360.), lim = c(-180.,360.)),
          axrat = list(lim = c(0.1,1.0), lim = c(0.1,1.0)),
          box = list(lim = c(-1.0,1.0), lim = c(-1.0,1.0))
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
    Data <- profitSetupData(image = image, mask = mask, sigma = sigma, segim = segim,
                            psf = psf, modellist = modellist, tofit = tofit, tolog = tolog,
                            priors = priors, intervals = intervals, magzero = ZP,
                            algo.func = 'optim', verbose = FALSE, like.func = like.func)

    ##########################
    ###### Initial Model #####
    ##########################
    
    png(file = paste(output_dir, 'initial_model.png', sep = "/"), width = 8, height = 2.5, units = "in", res = 300)
      profitLikeModel(parm = Data$init, Data = Data, makeplots = TRUE, whichcomponents = list(sersic = "all"))
    dev.off()
    
    png(file = paste(output_dir,'initial_1D_profile.png',sep = "/"), width = 7, height = 5, units = "in", res = 300)
      SBprof = GRAFitEllipsePlot(Data = Data, modellist = Data$modellist, pixscale = pix_scale, FWHM = 0.09, SBlim = 26, GRAFitlib = GRAFitlib)
    dev.off()
    

  } else if (nComp == 1) {
    ang  <- main_src$ang
    xcen <- main_src$xcen
    ycen <- main_src$ycen
    # mag <- -2.5*log10(c(main_src$flux*0.1,main_src$flux*0.9))-48.6
    # mag <- main_src$mag
    # flux <- main_src$flux
    mag  <- -2.5*log10(main_src$flux) + ZP
    re=sqrt(main_src$N50/(pi*main_src$axrat))  # = R50/pix_scale
    axrat <- main_src$axrat
    
    modellist = list(
      sersic = list(
        xcen = xcen,
        ycen = ycen,
        mag = mag ,
        re = re,
        nser = 1/main_src$con,
        ang = ang, 
        axrat = axrat, #min/maj: 1= o, 0= |
        box = 0     
      )
    )
    
    tofit = list(
      sersic = list(
        xcen = TRUE, # fit for xcen of the single component.
        ycen = TRUE, # fit for ycen of the single component.
        mag = TRUE, # Fit for both
        re = TRUE, # Fit for both
        nser = TRUE, # Fit for bulge
        ang = TRUE, # Fit for disk
        axrat = TRUE, # Fit for disk
        box = FALSE # Fit for neither
      )
    )
    
    # What parameters should be fitted in log space:
    
    tolog = list(
      sersic = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = FALSE,
        re = TRUE,    #re is best fit in log space
        nser = TRUE,  #nser is best fit in log space
        ang = FALSE,
        axrat = TRUE, #axrat is best fit in log space
        box = FALSE
      )
    )
    
    sigmas = c(2,2,5,1,1,30,0.3,Inf)
    
    sigmas = list(
      sersic = list(
        xcen = sigmas[1],
        ycen = sigmas[2],
        mag = sigmas[3],
        re = sigmas[4],
        nser = sigmas[5],
        ang = sigmas[6],
        axrat = sigmas[7],
        box = sigmas[8]
      )
    )
    
    priors = profitMakePriors(modellist, sigmas, tolog, allowflat = TRUE)
    
    # The hard intervals should also be specified in log space if relevant:
    
    intervals = list(
      sersic = list(
        xcen = list(lim = c(xcen-10, xcen+10)),
        ycen = list(lim = c(ycen-10, ycen+10)),
        mag = list(lim = c(10., 40.)),
        re = list(lim = c(1.0,dim(image)[1])),
        nser = list(lim = c(0.5, 15)),
        ang = list(lim = c(-180., 360.)),
        axrat = list(lim = c(0.1, 1.)),
        box = list(lim = c(-1, 1.))
      )
    )
    
    #tempCL=profitOpenCLEnv()
    Data <- profitSetupData(image = image, mask = mask, sigma = sigma, segim = segim,
                            psf = psf, modellist = modellist, tofit = tofit, tolog = tolog,
                            priors = priors, intervals = intervals, magzero = ZP,
                            algo.func = 'optim', verbose = FALSE, like.func = like.func)
    
    ##########################
    ###### Initial Model #####
    ##########################
    
    png(file = paste(output_dir, 'initial_model.png', sep = "/"), width = 8, height = 2.5, units = "in", res = 300)
      profitLikeModel(parm = Data$init, Data = Data, makeplots = TRUE, whichcomponents = list(sersic = 1))
    dev.off()
    
    png(file = paste(output_dir,'initial_1D_profile.png',sep = "/"),width=7,height=5,units="in",res=300)
      SBprof = GRAFitEllipsePlot(Data = Data, modellist = GRAFitAddFakeBulge(model =  Data$modellist, zeropoint = ZP),
                                pixscale = pix_scale, FWHM = 0.09, SBlim = 26, GRAFitlib = GRAFitlib) 
    dev.off()
  } # end of if on fit components.

  return(list(Data = Data, SBprof = SBprof))
  
}

# END


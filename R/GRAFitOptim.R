GRAFitOptim <- function( output_dir = output_dir, GRAFitlib = GRAFitlib, Model = profitLikeModel, nComp = 2,
                      Initial.Values = NULL, Data = Data, Algorithm = 'BFGS', zeropoint = ZP,
                      verbose = TRUE, main_src = main_src, ... ) {

  if( verbose ) cat("Running optimfit ......", '\n')

  source(paste(GRAFitlib,'/GRAFitAddFakeBulge.R', sep=''))
  source(paste(GRAFitlib,'/GRAFitEllipsePlot.R',sep=''))
  
  optimfit = optim(par = Data$init, Model, method = Algorithm,
                   Data = Data, control = list(fnscale = -1)) #,parscale=sigmas[which(unlist(tofit))]   method='L-BFGS-B'

  modeloptim = profitRemakeModellist(parm = optimfit$par, Data = Data)$modellist

  png(file = paste( output_dir, 'optim.png', sep = "/" ), width = 8, height = 2.5, units = "in", res = 200 )
    profitLikeModel( optimfit$par, Data, makeplots = TRUE, whichcomponents = list(sersic = 'all') )
  dev.off()

  png(file = paste( output_dir, 'optim_chisq.png', sep = "/"), width = 10, height = 7, units = "in", res = 200)
    profitLikeModel(optimfit$par, Data, makeplots = TRUE, whichcomponents = list(sersic = 'all'), plotchisq = T)
  dev.off()

    png(file = paste(output_dir, 'optim_1D_profile.png', sep = "/"), width = 10, height = 7, units = "in", res = 200)
      if (nComp == 1){
        fakemodellist = GRAFitAddFakeBulge( model =  modeloptim, zeropoint = zeropoint )
        SBprof = GRAFitEllipsePlot( Data = Data, modellist = fakemodellist, GRAFitlib = GRAFitlib, ...)
      } else {
        SBprof = GRAFitEllipsePlot( Data = Data, modellist = modeloptim, GRAFitlib = GRAFitlib, ...)
      }
    dev.off()

  #   if (nComp == 1){
  #   fakemodellist = GRAFitAddFakeBulge( model =  modeloptim, zeropoint = ZP )
  #   try(profitEllipsePlot(Data = Data, modellist = fakemodellist, 
  #                         pixscale = pix_scale, FWHM = FWHM, SBlim = SBlim))
  # } else if (nComp == 2){
  #   try(profitEllipsePlot(Data = Data, modellist = modeloptim, 
  #                         pixscale = pix_scale, FWHM = FWHM, SBlim = SBlim))
  #   
  # }
  #   dev.off()
  

  if(verbose) cat("optimfit finished ......", '\n')
  log8 = " optimfit:: DONE :D "

  return( list(optimfit = optimfit, modeloptim = modeloptim, SBprof = SBprof) )

}


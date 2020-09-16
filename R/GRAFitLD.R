GRAFitLD <- function( output_dir = output_dir, GRAFitlib = GRAFitlib, Model = profitLikeModel, nComp = 2,
                      Initial.Values = NULL, Data = Data, iteration = 1e4, zeropoint = ZP, pixscale = 0.03, 
                      FWHM = 0.09, SBlim = 26, Algorithm = 'CHARM', verbose = TRUE, ... ) {

  source(paste(GRAFitlib,'/GRAFitAddFakeBulge.R', sep=''))
  source(paste(GRAFitlib,'/GRAFitEllipsePlot.R',sep=''))
  source(paste(GRAFitlib,'/GRAFitri.R',sep=''))
  
  if( verbose ) cat("*** Performing FULL MCMC: Laplaces Demon ......", '\n')

  Data$algo.func="LD"   # if don't change the algorithm you'll get the ERROR: "Model must return a list.true"
  if ( is.null(Initial.Values) ) {Initial.Values = Data$init}
  # LDfit = LaplacesDemon(profitLikeModel, Initial.Values = Initial.Values, Data = Data, Status = 1000,
  #                       Iterations = iteration, Algorithm = Algorithm, Thinning = 1, Specs = list(alpha.star = 0.44))  #c(rep(0,14))

  LDfit = convergeFit(Data)      # from AllStarFit package by Dan Taranu.
  
  png(file = paste(output_dir,'LD_post.png',sep = "/"),width=30,height=30,units="in", res = 80)
    BestLD = GRAFitri(LDfit$Posterior1, samples = 1000, samptype = 'end')
  dev.off()

  png(file = paste(output_dir,'LD.png',sep = "/"), width = 8, height = 2.5,units = "in",res = 200)
    profitLikeModel( BestLD[,1], Data, makeplots = TRUE, whichcomponents = list(sersic = 'all'))
  dev.off()

  png(file = paste(output_dir,'LD_chisq.png',sep = "/"), width = 10, height = 7, units="in", res = 200)
    profitLikeModel( BestLD[,1], Data, makeplots = TRUE, whichcomponents = list(sersic='all'), plotchisq = T)
  dev.off()

  modelLD = profitRemakeModellist( BestLD[,1], Data$modellist, Data$tofit, Data$tolog )$modellist

    png(file = paste(output_dir,'/LD_1D_profile.png', sep = ''), width = 10, height = 7, units = "in", res = 200)
      if (nComp == 1){
        fakemodellist = GRAFitAddFakeBulge( model =  modelLD, zeropoint = zeropoint )
        SBprof = GRAFitEllipsePlot( Data = Data, modellist = fakemodellist, 
                                     pixscale = pixscale, FWHM = FWHM, 
                                     SBlim = SBlim, GRAFitlib = GRAFitlib)
      } else {
        SBprof = GRAFitEllipsePlot( Data = Data, modellist = modelLD, 
                                     pixscale = pixscale, FWHM = FWHM, 
                                     SBlim = SBlim, GRAFitlib = GRAFitlib)
      }
    dev.off()
    
  # png(file = paste(output_dir,'LD_post1.png',sep = "/"), width = 12, height = 12, units = "in", res = 300)
  #   magtri(LDfit$Posterior1[,3:ncol(LDfit$Posterior1)], samples = 1000, samptype =  'end')
  # dev.off()


  return(list( LDfit = LDfit, modelLD = modelLD, SBprof = SBprof ))

}
 # END

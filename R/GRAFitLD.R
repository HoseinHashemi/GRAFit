#' GRAFit: perform MCMC using LaplacesDemon
#'
#' @description This low level function acts as a wrapper around the \code{LalacesDemon} package to perform an MCMC sampling.
#' @param output_dir The directory to save outputs including the plots of the final model and standard \code{ProFit} plots.
#' @param Model profitLikeModel
#' @param nComp Number of components to be used in fitting. Default = 2 for a bulge+disk model. Alternatively could be 1 for a single Sersic model.
#' @param Initial.Values A vector of initial values equal in length to the number of parameters.
#' @param Data Data of class profit.data. The standard Data structure containing all inputs, images, modelsetup etc.
#' @param iteration The iterations for \code{MCMC} optimization if \code{optimMode = 'MCMC'}. Default = 10000
#' @param zeropoint Numeric scalar; the magnitude zero point.
#' @param pixscale Pixel scale in arcsecond/pixel.
#' @param FWHM The full width half max of the PSF in units of arc seconds. A vertical line is drawn at half this number (since we are plotting radius). The fits inside of the region is inherently hard because it is well within the PSF convolution kernel
#' @param SBlim 5 sigma surface brightness limit of the data typically in units of mag/asec^2. Default = 26.
#' @param Algorithm The MCMC algorithm to be used for the optimization. Default = 'CHARM'. You probably don't need to change this as CHARM is shown to be very robust algorithm.
#' @param verbose verbose.
#' @return Plots and optimized model.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitLA}}, \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @examples -
#' @export

GRAFitLD <- function( output_dir = output_dir, GRAFitlib = GRAFitlib, Model = profitLikeModel, nComp = 2,
                      Initial.Values = NULL, Data = Data, iteration = 1e4, zeropoint = ZP, pixscale = 0.03,
                      FWHM = 0.09, SBlim = 26, Algorithm = 'CHARM', verbose = TRUE, ... ) {

  if( verbose ) cat("*** Performing FULL MCMC: Laplaces Demon ......", '\n')

  Data$algo.func="LD"   # if don't change the algorithm you'll get the ERROR: "Model must return a list.true"
  if ( is.null(Initial.Values) ) {Initial.Values = Data$init}
  LDfit = LaplacesDemon(profitLikeModel, Initial.Values = Initial.Values, Data = Data, Status = 1000,
                        Iterations = iteration, Algorithm = Algorithm, Thinning = 1, Specs = list(alpha.star = 0.44))  #c(rep(0,14))

  # LDfit = AllStarFit::convergeFit(Data)      # from AllStarFit package by Dan Taranu.

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
        fakemodellist = GRAFitAddFakeBulge( model =  modelLD,
                                            zeropoint = zeropoint )

        SBprof = GRAFitEllipsePlot( Data = Data,
                                    modellist = fakemodellist,
                                    pixscale = pixscale,
                                    FWHM = FWHM,
                                    SBlim = SBlim,
                                    GRAFitlib = GRAFitlib)
      } else {
        SBprof = GRAFitEllipsePlot( Data = Data,
                                    modellist = modelLD,
                                    pixscale = pixscale,
                                    FWHM = FWHM,
                                    SBlim = SBlim,
                                    GRAFitlib = GRAFitlib)
      }
    dev.off()

  # png(file = paste(output_dir,'LD_post1.png',sep = "/"), width = 12, height = 12, units = "in", res = 300)
  #   magtri(LDfit$Posterior1[,3:ncol(LDfit$Posterior1)], samples = 1000, samptype =  'end')
  # dev.off()


  return(list( LDfit = LDfit, modelLD = modelLD, SBprof = SBprof ))

}
 # END

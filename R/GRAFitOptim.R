#' GRAFit: perform optim fit
#'
#' @description This low level function acts as a wrapper around the \code{\link[stats]{optim}}. It can handle a few optimization algorithms including: \code{"Nelder-Mead"}, \code{"CG"}, \code{"BFGS"}, \code{"L-BFGS-B"}, \code{"SANN"}, \code{"Brent"}.
#' @param output_dir The directory to save outputs including the plots of the final model and standard \code{ProFit} plots.
#' @param Model profitLikeModel
#' @param nComp Number of components to be used in fitting. Default = 2 for a bulge+disk model. Alternatively could be 1 for a single Sersic model.
#' @param Initial.Values Initial values.
#' @param Data Data of class profit.data. The standard Data structure containing all inputs, images, modelsetup etc.
#' @param Algorithm A string that specifies the algorithm for optimization. Default = \code{'BFGS'} a quasi-Newton method introduced by Broyden, Fletcher, Goldfarb and Shanno. Other options \code{"Nelder-Mead"}, \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"}, \code{"Brent"}. See \code{\link[stats]{optim}} for more details.
#' @param zeropoint Numeric scalar; the magnitude zero point.
#' @param verbose verbose.
#' @param ... extra parameters to be parsed to \code{\link[GRAFit]{GRAFitEllipsePlot}}.
#' @return A list of optimized fit containing both the output of \code{LaplaceApproximation} and \code{profitRemakeModellist}.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitLD}}, \code{\link[GRAFit]{GRAFitLA}}, \code{\link[GRAFit]{GRAFitcutoutWCS}}
#' @examples -
#' @export



GRAFitOptim <- function( output_dir = output_dir, GRAFitlib = GRAFitlib, Model = profitLikeModel, nComp = 2,
                      Initial.Values = NULL, Data = Data, Algorithm = 'BFGS', zeropoint = ZP,
                      verbose = TRUE, ... ) {

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

  if(verbose) cat("optimfit finished ......", '\n')
  log8 = " optimfit:: DONE :D "

  return( list(optimfit = optimfit, modeloptim = modeloptim, SBprof = SBprof) )

}


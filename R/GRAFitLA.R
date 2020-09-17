#' GRAFit: perform Laplace Approximation
#'
#' @description This low level function acts as a wrapper around the \code{\link[LaplacesDemon]{LaplaceApproximation}}.
#' @param Model profitLikeModel
#' @param Initial.Values Initial values.
#' @param Data Data of class profit.data. The standard Data structure containing all inputs, images, modelsetup etc.
#' @param iteration An integer. The iterations for optimization. Default = 1000
#' @param Method A string that specifies the mothod for Laplace Approximation. Default = \code{'LM'} for Levenberg-Marquardt. See \code{\link[LaplacesDemon]{LaplaceApproximation}} for all options.
#' @param verbose verbose.
#' @return A list of optimized fit containing both the output of \code{LaplaceApproximation} and \code{profitRemakeModellist}.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitLD}}, \code{\link[LaplacesDemon]{LaplaceApproximation}}
#' @examples -
#' @export

GRAFitLA <- function( Model = profitLikeModel, Initial.Values = Data$init, Data = Data,
                      iteration = 1e3, Method = 'LM', verbose = TRUE ) {

  if ( verbose ) cat("Running Laplace Approximation ......", '\n')
  #       pdf(file = paste(output_dir,'/LAfit.pdf',sep = ""), width = 10, height = 8)
  Data$algo.func = "LA"   # if don't change the algorithm you'll get the ERROR: "Model must ereturn a list.true"

  LAfit = LaplaceApproximation( Model =  Model, parm = Initial.Values, Data = Data, Iterations = iteration,
                             Method = Method, CovEst = 'Identity', sir  = FALSE)
  #       profitLikeModel(LAfit$Summary1[,1],Data,makeplots=TRUE,whichcomponents=list(sersic='all'))
  #       profitLikeModel(LAfit$Summary1[,1],Data,makeplots=TRUE,whichcomponents=list(sersic='all'),plotchisq = T)
  modelLA = profitRemakeModellist( LAfit$Summary1[,1], Data$modellist, Data$tofit, Data$tolog)$modellist
  #       if (nComps == 1){
  #         try(profitEllipsePlot(Data=Data,modellist=GRAFit_add_fake_bulge(modelLA),pixscale=pix_scale,FWHM=0.124,SBlim=26))
  #       } else if (nComps == 2){
  #         try(profitEllipsePlot(Data=Data,modellist=modelLA,pixscale=pix_scale,FWHM=0.124,SBlim=26))
  #       }
  #       dev.off()
  log9=" Laplace Approximation:: DONE :D "

  return(list( LAfit = LAfit, modelLA = modelLA ))
}


#' GRAFit: Find the Target Object in the Cutout.
#'
#' @description This function finds the main target/galaxy a cutout.
#' @param src_list The output of \code{\link[ProFound]{ProFound}} source detection.
#' @param imDim An n*m matrix. The dimension of the cutout.
#' @param main_pin Should the location of the main object on the cutout be shown by a symbol. This can only be \code{TRUE} if the image is open.
#' @return A file containing the coordination of each frame's center.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitDynamo_v2}}
#' @examples -
#' @export

GRAFitMainFinder = function(src_list = NULL,
                            imDim = NULL,
                            main_pin = FALSE) {

  nSrc = length(src_list$segID) # number of detected sources.

  x0 = imDim[1]/2; y0 = imDim[2]/2 # position of cut-out centre.
  sepArr = array(0, dim=nSrc) # empty separation array

  # Calculate separation of sources to the centres of image.
  for (i in seq(1,nSrc)) {
    sepArr[i] = sqrt( (src_list$xcen[i] - x0)^2 + (src_list$ycen[i] - y0)^2 )
  }

  mainObj = which.min(sepArr) # Main source is the one with the minimum separation from the cutout centre
  main = src_list[mainObj, ]

  if (main_pin) points(main$xmax, main$ymax, pch= 4, col = 'red')

  return( main )
}

# END

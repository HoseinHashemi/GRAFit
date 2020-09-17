#' GRAFit: Overplot the Ellipsoidal Aperture.
#'
#' @description This function overplots the ellipsoidal apertures on the image.
#' @param x,y The location of the aperture to be plotted.
#' @param a,b minor and major axes of the apperture.
#' @param angle Position angle
#' @param border The colour of the ellipse.
#' @return Overplots the ellipsoidal apertur on the image.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitEllipsePlot}}
#' @examples -
#' @export
#'


GRAFitAperturePlot <- function(x = x, y = y, a = a, b = b, angle = angle, border = 'red') {
  points(x, y, pch = 4, col = border)
  draw.ellipse(x,y,a,b,angle, border = border)
}

# END

#' GRAFit: calculates the central surface brightness
#'
#' @description This function calculates the central surface brightness of a galaxy from its Sersic index profile using the euation in Graham & Driver 2005.
#' @param n Sersic index
#' @param mag Magnitude
#' @param Re The effective radius either in pixel for which the PixScale should be provided or arcsecond.
#' @param PixScale Pixel scale in arcsecond/pixel. \code{PixScale = 1} (default) means that Re is in the units of arcsecond. Should be necessarily provided if the Re is in units of pixels.
#' @return The central surface brightness of galaxy in units of mag/asec^2
#' @author Hosein Hashemizadeh
#' @references Graham A. W., Driver S. P., 2005, PASA, 22, 118
#' @seealso \code{\link[GRAFit]{GRAFitEllipsePlot}}
#' @examples -
#' GRAFitMu0(n = 2, mag = 21, Re = 20, PixScale = 0.03)
#' @export


GRAFitMu0 <- function(n, mag, Re, PixScale = 1) {

  Re = Re * PixScale
  bn = 1.9992 * n - 0.3271

  mu0 = mag + 5 * log10(Re) - 2.5 * log10( bn ** (2*n) / (pi * gamma(2 * n + 1))  )

  return(mu0)
}

#END

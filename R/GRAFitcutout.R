#' GRAFit: Image CutoutWCS
#'
#' @description This function cuts a subset of an image with associated WCS system. This is originally developed as part of \code{\link[magicaxis]{magicaxis}} package. It has been modified to handle different paddings with arbitrary values.
#' @param image	Numeric matrix; required, the image we want to decorate. If image is a list as created by readFITS, read.fits of magcutoutWCS then the image part of these lists is parsed to image and the correct header part is parsed to header.
#' @param header Full FITS header in table or vector format. Legal table format headers are provided by the read.fitshdr function or the hdr list output of read.fits in the astro package; the hdr output of readFITS in the FITSio package or the header output of magcutoutWCS. If a header is provided then key words will be taken from here as a priority. Missing header keywords are printed out and other header option arguments are used in these cases.
#' @param loc	Numeric vector; the [x,y] (magcutout) or [x,y] / [RA,Dec] (magcutoutWCS) location where we want to cutout the image. For magcutoutWCS the unit type can be specified with the loc.type option. Either it is WCS in degrees [RA,Dec] (coord) or pixel [x,y] of the image (image).
#' @param box	Numeric vector; the dimensions of the box to cut out from image centred on loc. For magcutout the box unit is always pixels. For magcutoutWCS the unit type can be specified with the loc.type option. Either it is pixels or asec (see loc.type).
#' @param shiftloc Logical; should the cutout centre shift from loc away from the image edge if the desired box size extends beyond the edge of the image?
#' @param paddim	Logical; should the cutout be padded with image data until it meets the desired box size (if shiftloc is TRUE) or padded with NAs for data outside the image boundary otherwise?
#' @param padVal the valuse to be replaced if a cutout needed to be padded.
#' @param plot Logical; should a plot be generated?
#' @param CRVAL1 FITS header CRVAL1 for the CTYPE1 projection system. This is the RA in degrees at the location of CRPIX1.
#' @param CRVAL2	FITS header CRVAL2 for the CTYPE2 projection system. This is the Dec in degrees at the location of CRPIX2.
#' @param CRPIX1	FITS header CRPIX1 for the CTYPE1 projection system. This is the x pixel value at the location of CRVAL1.
#' @param CRPIX2	FITS header CRPIX2 for the CTYPE2 projection system. This is the y pixel value at the location of CRVAL2.
#' @param CD1_1	FITS header CD1_1 for the CTYPE1 projection system. Change in CTYPE1 in degrees along x-Axis.
#' @param CD1_2	FITS header CD1_2 for the CTYPE1 projection system. Change in CTYPE1 in degrees along y-Axis.
#' @param CD2_1	FITS header CD2_1 for the CTYPE2 projection system. Change in CTYPE2 in degrees along x-Axis.
#' @param CD2_2	FITS header CD2_2 for the CTYPE2 projection system. Change in CTYPE2 in degrees along y-Axis.
#' @param CTYPE1	The RA projection system type. Either 'RA–TAN' for Tan Gnomonic (default), or 'RA–SIN' for Sine Orthographic. 'RA–NCP' is approximated by Sine Orthographic with a warning. Over-ridden by the FITS header.
#' @param CTYPE2	The DEC projection system type. Either 'DEC–TAN' for Tan Gnomonic (default), or 'DEC–SIN' for Sine Orthographic. 'DEC–NCP' is approximated by Sine Orthographic with a warning. Over-ridden by the FITS header.
#' @param coord.type	The units of loc for magcutoutWCS. Allowed options are 'deg' for degress and 'sex' for sexigesimal (i.e. HMS for RA and DMS for Deg).
#' @param sep	When coord.type='sex', sep defines the type of separator used for the HMS and DMS strings (i.e. H:M:S and D:M:S would be sep=':', which is the default). See hms2deg and dms2deg for more details.
#' @param loc.type	Character vector; specifies what type of location is being provided. The first element specifies the coordinate type for loc and the second element is the coordinate type for box. Either it is WCS in degrees [RA,Dec] / asec ('coord') or pixel [x,y] of the image ('image'). If only one element is provided then the same coordinate type is used for both loc and box.
#' @param ... Dots are parsed to either \code{\link[magimage]{magimage}} (GRAFitcutout) or \code{\link[magimage]{magimageWCS}} (GRAFitcutoutWCS)
#' @return A file containing the coordination of each frame's center.
#' @author Aaron Robotham; Modified by: Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitcutout}}
#' @examples -
#' @export
#'
GRAFitcutout <- function (image, loc = dim(image)/2, box = c(100, 100), shiftloc = FALSE,
          paddim = TRUE, padVal = 0, plot = FALSE, ...)
{
  loc = as.numeric(loc)
  xcen = loc[1]
  ycen = loc[2]
  loc = ceiling(loc)
  xlo = ceiling(loc[1] - (box[1]/2 - 0.5))
  xhi = ceiling(loc[1] + (box[1]/2 - 0.5))
  ylo = ceiling(loc[2] - (box[2]/2 - 0.5))
  yhi = ceiling(loc[2] + (box[2]/2 - 0.5))
  loc.diff = c(x = xlo - 1, y = ylo - 1)
  expand = paddim && shiftloc
  diffxlo = xlo - 1
  if (diffxlo < 0) {
    xlo = 1
    if (expand)
      xhi = xlo + (box[1] - 1)
  }
  diffxhi = xhi - dim(image)[1]
  if (diffxhi > 0) {
    xhi = dim(image)[1]
    if (expand) {
      xlo = xlo - diffxhi
      if (xlo < 1)
        xlo = 1
    }
  }
  diffylo = ylo - 1
  if (diffylo < 0) {
    ylo = 1
    if (expand)
      yhi = ylo + (box[2] - 1)
  }
  diffyhi = yhi - dim(image)[2]
  if (diffyhi > 0) {
    yhi = dim(image)[2]
    if (expand) {
      ylo = ylo - diffyhi
      if (ylo < 1)
        ylo = 1
    }
  }
  if (!paddim && !shiftloc) {
    if (diffxlo < 0 && (-diffxlo > diffxhi))
      xhi = xhi - max(diffxhi, 0) + diffxlo
    if (diffxhi > 0 && (-diffxlo < diffxhi))
      xlo = xlo + diffxhi - min(diffxlo, 0)
    if (diffylo < 0 && (-diffylo > diffyhi))
      yhi = yhi - max(diffyhi, 0) + diffylo
    if (diffyhi > 0 && (-diffylo < diffyhi))
      ylo = ylo + diffyhi - min(diffylo, 0)
  }
  xsel = xlo:xhi
  ysel = ylo:yhi
  if (xsel[2] == 0 | ysel[2] == 0) {
    image = matrix(padVal, box[1], box[2])
  }
  else {
    image = image[xsel, ysel]
    if (paddim && !shiftloc && any(c(diffxlo, -diffxhi, diffylo,
                                     -diffyhi) < 0)) {
      padded = matrix(padVal, box[1], box[2])
      padded[xsel - diffxlo, ysel - diffylo] = image
      image = padded
    }
  }
  if (shiftloc) {
    loc.diff = c(x = xlo - 1, y = ylo - 1)
  }
  output = list(image = image, loc = c(x = xcen - xlo + 1,
                                       y = ycen - ylo + 1), loc.orig = c(x = xcen, y = ycen),
                loc.diff = loc.diff, xsel = xsel, ysel = ysel)
  if (plot) {
    if (all(is.na(image))) {
      image[] = 0
      magimage(image, ...)
    }
    else {
      magimage(image, ...)
    }
  }
  return = output
}

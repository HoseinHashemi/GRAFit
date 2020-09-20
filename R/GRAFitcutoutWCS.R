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

GRAFitcutoutWCS <- function ( image,
                              header,
                              loc,
                              box = c(100, 100),
                              shiftloc = FALSE,
                              paddim = TRUE,
                              padVal = 0,
                              plot = FALSE,
                              CRVAL1 = 0,
                              CRVAL2 = 0,
                              CRPIX1 = 0,
                              CRPIX2 = 0,
                              CD1_1 = 1,
                              CD1_2 = 0,
                              CD2_1 = 0,
                              CD2_2 = 1,
                              coord.type = "deg",
                              sep = ":",
                              loc.type = c("coord", "coord"), ...)
{

  if (length(loc.type) == 1) {
    loc.type = rep(loc.type, 2)
  }
  if (!missing(image)) {
    if (any(names(image) == "imDat") & missing(header)) {
      imtype = "FITSio"
      header = image$hdr
      image = image$imDat
    }
    if (any(names(image) == "dat") & missing(header)) {
      imtype = "astro"
      header = image$hdr[[1]]
      header = data.frame(key = header[, 1], value = header[,
                                                            2], stringsAsFactors = FALSE)
      image = image$dat[[1]]
    }
    if (any(names(image) == "image") & missing(header)) {
      header = image$header
      image = image$image
      if (is.matrix(header) | is.data.frame(header)) {
        imtype = "astro"
      }
      else {
        imtype = "FITSio"
      }
    }
    if (!missing(header)) {
      if (is.matrix(header) | is.data.frame(header)) {
        imtype = "astro"
      }
      else {
        imtype = "FITSio"
      }
    }
  }
  if (missing(loc)) {
    loc = magWCSxy2radec(dim(image)[1]/2, dim(image)[2]/2,
                         header = header, CRVAL1 = CRVAL1, CRVAL2 = CRVAL2,
                         CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, CD1_1 = CD1_1,
                         CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2)[1, ]
    tempxy = cbind(dim(image)[1]/2, dim(image)[2]/2)
  }
  else {
    if (loc.type[1] == "coord") {
      if (coord.type == "sex") {
        loc[1] = hms2deg(loc[1], sep = sep)
        loc[2] = dms2deg(loc[2], sep = sep)
      }
      loc = as.numeric(loc)
      tempxy = magWCSradec2xy(loc[1], loc[2], header = header,
                              CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1,
                              CRPIX2 = CRPIX2, CD1_1 = CD1_1, CD1_2 = CD1_2,
                              CD2_1 = CD2_1, CD2_2 = CD2_2)
    }
    else if (loc.type[1] == "image") {
      tempxy = rbind(loc)
      loc = magWCSxy2radec(loc[1], loc[2], header = header,
                           CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1,
                           CRPIX2 = CRPIX2, CD1_1 = CD1_1, CD1_2 = CD1_2,
                           CD2_1 = CD2_1, CD2_2 = CD2_2)[1, ]
    }
  }
  xcen = tempxy[1, 1]
  ycen = tempxy[1, 2]
  if (loc.type[2] == "coord") {
    box = box/3600
    tempxy = magWCSradec2xy(loc[1] - box[1]/2/cos(loc[2] *
                                                    pi/180), loc[2], header = header, CRVAL1 = CRVAL1,
                            CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2,
                            CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2)
    xlo = xcen - sqrt((tempxy[1, 1] - xcen)^2 + (tempxy[1,
                                                        2] - ycen)^2)
    tempxy = magWCSradec2xy(loc[1] + box[1]/2/cos(loc[2] *
                                                    pi/180), loc[2], header = header, CRVAL1 = CRVAL1,
                            CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2,
                            CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2)
    xhi = xcen + sqrt((tempxy[1, 1] - xcen)^2 + (tempxy[1,
                                                        2] - ycen)^2)
    tempxy = radec2xy(loc[1], loc[2] - box[2]/2, header = header,
                      CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1,
                      CRPIX2 = CRPIX2, CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1,
                      CD2_2 = CD2_2)
    ylo = ycen - sqrt((tempxy[1, 1] - xcen)^2 + (tempxy[1,
                                                        2] - ycen)^2)
    tempxy = magWCSradec2xy(loc[1], loc[2] + box[2]/2, header = header,
                            CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1,
                            CRPIX2 = CRPIX2, CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1,
                            CD2_2 = CD2_2)
    yhi = ycen + sqrt((tempxy[1, 1] - xcen)^2 + (tempxy[1,
                                                        2] - ycen)^2)
    xtemp = sort(c(xlo, xhi))
    xlo = ceiling(xtemp[1])
    xhi = ceiling(xtemp[2])
    ytemp = sort(c(ylo, yhi))
    ylo = ceiling(ytemp[1])
    yhi = ceiling(ytemp[2])
    box = c(xhi - xlo + 1, yhi - ylo + 1)
  }
  else {
  }
  cutout = GRAFitcutout(image, loc = c(xcen, ycen), box = box,
                     shiftloc = shiftloc, paddim = paddim, padVal = padVal, plot = FALSE)
  cut_image = cutout$image
  xlo = cutout$loc.diff[1] + 1
  xhi = xlo + dim(cut_image)[1] - 1
  ylo = cutout$loc.diff[2] + 1
  yhi = ylo + dim(cut_image)[2] - 1
  xcen.new = xcen - xlo + 1
  ycen.new = ycen - ylo + 1
  pixscale = getpixscale(header = header, CD1_1 = CD1_1, CD1_2 = CD1_2,
                         CD2_1 = CD2_1, CD2_2 = CD2_2)
  loc.diff = c(xlo - 1, ylo - 1)
  cut_xlo = 1
  cut_xhi = dim(cut_image)[1]
  cut_ylo = 1
  cut_yhi = dim(cut_image)[2]
  usr.WCS = rbind(magWCSxy2radec(xlo - 1, ylo - 1, header = header,
                                 CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2,
                                 CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2),
                  magWCSxy2radec(xlo - 1, yhi, header = header, CRVAL1 = CRVAL1,
                                 CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2,
                                 CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2),
                  magWCSxy2radec(xhi, ylo - 1, header = header, CRVAL1 = CRVAL1,
                                 CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2,
                                 CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2),
                  magWCSxy2radec(xhi, yhi, header = header, CRVAL1 = CRVAL1,
                                 CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2,
                                 CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2))
  usr.WCS = cbind(x.cut = c(cut_xlo - 1, cut_xlo - 1, cut_xhi,
                            cut_xhi), y.cut = c(cut_ylo - 1, cut_yhi, cut_ylo - 1,
                                                cut_yhi), x.orig = c(xlo - 1, xlo - 1, xhi, xhi), y.orig = c(ylo -
                                                                                                               1, yhi, ylo - 1, yhi), usr.WCS)
  approx.map.RA = approxfun(seq(usr.WCS[1, "RA"], usr.WCS[4,
                                                          "RA"], len = 100), seq(usr.WCS[1, "x.cut"], usr.WCS[4,
                                                                                                              "x.cut"], len = 100))
  approx.map.Dec = approxfun(seq(usr.WCS[1, "Dec"], usr.WCS[4,
                                                            "Dec"], len = 100), seq(usr.WCS[1, "y.cut"], usr.WCS[4,
                                                                                                                 "y.cut"], len = 100))
  approx.map = function(RA, Dec) {
    if (length(dim(RA)) == 2) {
      Dec = RA[, 2]
      RA = RA[, 1]
    }
    return = cbind(x = approx.map.RA(RA), y = approx.map.Dec(Dec))
  }
  if (!missing(header)) {
    dimdiff = dim(cut_image) - dim(image)
    hdradd = list(CRPIX1 = -loc.diff[1], CRPIX2 = -loc.diff[2],
                  NAXIS1 = dimdiff[1], NAXIS2 = dimdiff[2])
    if (imtype == "FITSio") {
      for (hdrname in names(hdradd)) {
        if (hdradd[[hdrname]] != 0) {
          hdrrow = which(header == hdrname) + 1
          header[hdrrow] = as.character(as.numeric(header[hdrrow]) +
                                          hdradd[[hdrname]])
        }
      }
    }
    else if (imtype == "astro") {
      for (hdrname in names(hdradd)) {
        if (hdradd[[hdrname]] != 0) {
          hdrrow = which(header[, "key"] == hdrname)
          header[hdrrow, "value"] = as.character(as.numeric(header[hdrrow,
                                                                   "value"]) + hdradd[[hdrname]])
        }
      }
    }
    else {
      header = NULL
    }
  }
  else {
    header = NULL
  }
  output = list(image = cut_image, loc = c(x = as.numeric(xcen.new),
                                           y = as.numeric(ycen.new)), loc.orig = c(x = as.numeric(xcen),
                                                                                   y = as.numeric(ycen)), loc.diff = c(as.numeric(loc.diff[1]),
                                                                                                                       as.numeric(loc.diff[2])), xsel = xlo:xhi, ysel = ylo:yhi,
                loc.WCS = loc, scale.WCS = pixscale, usr.WCS = usr.WCS,
                approx.map = approx.map, header = header)
  if (plot) {
    if (all(is.na(cut_image))) {
      cut_image[] = 0
      magimageWCS(image = cut_image, header = header, CRVAL1 = CRVAL1,
                  CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2,
                  CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1,
                  CD2_2 = CD2_2, ...)
      cut_image[] = NA
    }
    else {
      magimageWCS(image = cut_image, header = header, CRVAL1 = CRVAL1,
                  CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2,
                  CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1,
                  CD2_2 = CD2_2, ...)
    }
  }
  return = output
}

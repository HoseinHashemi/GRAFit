# Author: Aaron Robotham as part of the magicaxis package (original: magcutoutWCS).
# Modified by: Hosein Hashemizadeh as part of GRAFit package.

GRAFitcutoutWCS <- function ( GRAFitlib_path = GRAFitlib, image, header, loc, box = c(100, 100), shiftloc = FALSE, 
          paddim = TRUE, padVal = 0, plot = FALSE, CRVAL1 = 0, CRVAL2 = 0, CRPIX1 = 0, 
          CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, coord.type = "deg", 
          sep = ":", loc.type = c("coord", "coord"), ...) 
{
  source(paste(GRAFitlib,'/GRAFitcutout.R',sep=''))
  
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
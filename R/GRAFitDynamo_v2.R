#' GRAFit: Dynamic Cutout
#'
#' @description This high level function utilizes the effective radius (Re) of the object to dynamically make a cutout around it.
#' @param im Image; a 2D image matrix.
#' @param loc Objects location; \code{c(RA,DEC)}
#' @param init_xrad Initial cutout size in x-direction in units of asec.
#' @param init_xrad Initial cutout size in x-direction in units of asec.
#' @param R90 The radius containing 90\% of the total flux. In units of asec.
#' @param finalCut Size of the final cutout as finalCut*R90. Default = 5, i.e., 5*R90.
#' @param finalBoxCar Integer vector; the dimensions of the box car filter to estimate the sky with. Default = 3, i.e. 3*R90. For convenience, if length 1 then both dimensions of box used internally are assumed to equal the specified box. i.e. 200 would be interpreted as c(200,200). Dependent default arguments (grid, boxadd and skypixmin) are updated sensibly.
#' @param header Full FITS header in table or vector format.
#' @param pix_scale Pixel scale in units of arcsecond/pixel.
#' @param verbose verbose.
#' @param magzero Numeric scalar; the magnitude zero point.
#' @param tolerance Numeric scalar; the minimum height of the object in the units of skyRMS between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbours, which is the highest. The range 1-5 offers decent results usually. Passed to \code{\link[ProFound]{profoundMakeSegim}}.
#' @param smooth Logical; should smoothing be done on the target image? Passed to \code{\link[ProFound]{profoundMakeSegim}}.
#' @param sigma Numeric scalar; standard deviation of the blur used when smooth=TRUE. Passed to \code{\link[ProFound]{profoundMakeSegim}}.
#' @param sky User provided estimate of the absolute sky level. If this is not provided then it will be computed internally using \code{\link[ProFound]{profoundMakeSkyGrid}}. Can be a scalar or a matrix matching the dimensions of image (allows values to vary per pixel). This will be subtracted off the image internally, so only provide this if the sky does need to be subtracted!
#' @param pixcut Integer scalar; the number of pixels required to identify an object.
#' @param skycut Numeric scalar; the lowest threshold to make on the image in units of the skyRMS.
#' @param boundstats Logical; if TRUE then various pixel boundary statistics are computed (\code{Nedge}, \code{Nsky}, \code{Nobject}, \code{Nborder}, \code{edge_frac}, \code{edge_excess} and \code{FlagBorder}). If FALSE these return NA instead (saving computation time).
#' @param rotstats Logical; if TRUE then the \code{asymm}, \code{flux_reflect} and \code{mag_reflect} are computed, else they are set to NA. This is because they are very expensive to compute compared to other photometric properties.
#' @param ImPlot Logical; should the image be plotted?
#' @param output_dir The directory to save outputs.
#' @return Plots as well we matrices of the image, segmentation, sky etc.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitLA}}, \code{\link[LaplacesDemon]{LaplacesDemon}}
#' @examples -
#' @export

GRAFitDynamo_v2 <- function( im = NULL,
                             loc = c(ra_deg, dec_deg),
                             init_xrad = NULL,
                             init_yrad = NULL,
                             R90 = NULL,
                             finalCut = 5,
                             finalBoxCar = 3,
                             header = NULL,
                             pix_scale = NULL,
                             magzero = NULL,
                             tolerance = 7,
                             smooth = TRUE,
                             sigma = 7,
                             sky = NULL,
                             pixcut = 3,
                             skycut = 1.1,
                             boundstats = TRUE,
                             rotstats = TRUE,
                             ImPlot = TRUE,
                             output_dir = NULL,
                             verbose = TRUE,
                             ... ) {

  if(verbose) cat( paste( "*** Doing Dynamic cut-out ", sep = "" ),'\n' )

  if (!is.null(init_xrad)) {
    xrad = init_xrad; yrad = init_yrad  # in arcsec
  } else {
    rad = round(15*R90)   # in arcsec
  }

  xy = magWCSradec2xy(RA = loc[1],
                      Dec = loc[2],
                      header = header)

  cut_im =  GRAFitcutout( GRAFitlib_path = GRAFitlib,
                          image = im,
                          loc = c(xy[1],xy[2]),
                          box = c( ceiling(rad/0.03), ceiling(rad/0.03) ),
                          paddim = TRUE,
                          header = header,
                          grid.lwd = 1,
                          plot = ImPlot,
                          padVal = 0 )

  # cut_im = GRAFitcutoutWCS( GRAFitlib_path = GRAFitlib, image = image, loc = loc, box = c( xrad, yrad ),
  #                              paddim = TRUE, header = header, grid.lwd = 1, plot = ImPlot, padVal = 0 )
  #

  image <- cut_im$image

  invisible(
  ifelse ( (image[1,] == 0) | ( image[,1] == 0 ),
           image <- floodFill( cut_im$image, c(1,1), col = NA ),
           (image = image) )
  )
  invisible(
  ifelse ( ( image[dim(image)[1], ] == 0 ) | (image[ ,dim(image)[2] ] == 0 ),
           image <- floodFill( cut_im$image, c(dim(image)[1], dim(image)[2]), col = NA ),
           (image = image) )
  )

  # mask = matrix(0, dim(image)[1], dim(image)[1]); mask[which(is.na(image))] = 1
  # mask = matrix(0, dim(image)[1], dim(image)[1]); mask[which(image == 0)] = 1
  box_car = round(5*R90/pix_scale)

  if((box_car %% 2) != 0) {
    box_car = box_car-1
  }

  print("----- profound")
  print(paste("box = ", box_car, "R90 = ", R90))

  seg <- profoundProFound( image = image,
                           tolerance = tolerance,
                           sigma = sigma,
                           sky = sky,      # mask = image==0,
                           smooth = smooth,
                           pixcut = pixcut,
                           skycut = skycut,
                           magzero = magzero,
                           pixscale = pix_scale,
                           box = c(box_car,box_car),
                           type = "bicubic",
                           header = header,
                           boundstats = boundstats,
                           groupstats = FALSE,
                           rotstats = rotstats,
                           plot = FALSE, ... )  #gain = gain_ADU

  segstats <- seg$segstats
  main_src0 = GRAFitMainFinder(src_list =  segstats,
                               imDim = dim(image))

  print("----- dilation")

  seg_dilate = profoundMakeSegimDilate(image = image,
                                       seg$segim, size = 51, #expand = main_src$segID,
                                       boundstats = TRUE,
                                       rotstats = TRUE,
                                       header = header,
                                       pixscale = pix_scale,
                                       magzero = magzero,
                                       plot = FALSE )

  log004 = " dilation:: DONE :D "
  segim_dil <- seg_dilate$segim
  mask = segim_dil == 0
  segstats <- seg_dilate$segstats
  main_src = GRAFitMainFinder(src_list =  segstats,
                              imDim = dim(image))

  i = finalCut
    box = c(round(i*R90/pix_scale),
            round(i*R90/pix_scale))

    cut_im_final_orig <- magcutoutWCS( image = im,
                                       header = header,
                                       loc = loc,
                                       box = c(round(i*R90), round(i*R90)) )

    loc = c(main_src$xcen, main_src$ycen)
    box = dim(cut_im_final_orig$image)

    segim <- GRAFitcutout( GRAFitlib_path = GRAFitlib, image = segim_dil,
                              box = box )

    sky <- magcutout( seg$sky, box = box )

    skyRMS <- magcutout( seg$skyRMS, box = box )

    objects <- magcutout( seg_dilate$objects, box = box )

    box_car = round(finalBoxCar*main_src$R90/pix_scale)
    seg2 <- profoundProFound( image = cut_im_final_orig$image,
                              segim = segim$image ,
                              tolerance = tolerance,
                              sigma = sigma,
                              smooth = smooth ,
                              pixcut = pixcut,
                              skycut = skycut,
                              magzero = magzero,
                              pixscale = pix_scale,
                              box = c(box_car, box_car),
                              header = header,
                              boundstats = boundstats,
                              groupstats = FALSE,
                              rotstats = rotstats,
                              plot = ImPlot,
                              lowmemory = F, ... )

    main_src2 = GRAFitMainFinder(src_list =  seg2$segstats,
                                 imDim = dim(cut_im_final_orig$image))

  if (!is.null(output_dir)) {
    png(filename = paste(output_dir,"/stamps.png", sep = ''),
             width = 1200, height = 1275)
      par(mfrow = c(2,2), oma = c(4, 4, 4, 4))  #, cex.lab = 3, cex.axis = 3
      profoundSegimPlot(image, seg$segim,
                        xlab = "x/pix", ylab = "y/pix",
                        cex.lab = 3,
                        cex.axis = 3,
                        mgp = c(2.5,1.5,0),
                        mtline = 5)
      title(main = "Initial cut out = 15R90 (Y UltraVISTA)", cex.main = 2.5)
      profoundSegimPlot(cut_im_final_orig$image, segim$image,
                        xlab = "x/pix", ylab = "",
                        cex.lab = 3,
                        cex.axis = 3,
                        mgp = c(2.5,1.5,0),
                        mtline = 5)
      title(main = paste("Final cut out = ",i, "R90 (Y UltraVISTA)", sep = ''), cex.main = 2.5)
      magimage(seg$sky, xlab = "x/pix", ylab = "y/pix",
               cex.lab = 3,
               cex.axis = 3,
               mgp = c(2.5,1.5,0),
               mtline = 5)
      title(main = "sky", cex.main = 2.5)
      profoundSkyEst(image = image, objects = seg$objects_redo, plot = TRUE, xlab = "sky pix Val",
                     cex.lab = 2,
                     cex.axis = 2,
                     mgp = c(2.5,1.5,0),
                     mtline = 5)
      title(main = "sky pix val distribution", cex.main = 2.5)
    dev.off()

    log3 = " Image cutout:: DONE :D "
    output = return(list( orig_image = cut_im_final_orig$image, sky_red_image = cut_im_final_orig$image-sky$image,
                          segim = segim$image, main_src = main_src2, sky = sky$image, skyRMS = skyRMS$image,
                          objects = objects$image))

  } else {
    log3 = " Image cutout:: DONE :D "
    output = return(list( orig_image = cut_im_final_orig$image, sky_red_image = cut_im_final_orig$image-sky$image,
                          segim = segim$image, main_src = main_src2, sky = sky$image, skyRMS = skyRMS$image,
                          objects = objects$image))

  }

}

# END


#' GRAFit: Setup Tiny Tim Commands for PSF Generation
#'
#' @description Low-level function to setup the requiered Tiny Tim commands for HST PSF generation.
#' @param wrk_dir Working directory, where the GRAFit outputs will be saved.
#' @param output_dir The directory to save the PSF.
#' @param name The name of output file
#' @param CCDCHIP Which ACS's CCD chip the object is located in. Allowed values: 1 and 2.
#' @param x,y target_loc unit; CCD coordination.
#' @param filter What filter is your imaging data. E.g., \code{"f814w"}. If missing it will be read from the \code{header}.
#' @param focus HST focus values.
#' @param ebmv B-V extinction.
#' @param jitter Telescope's jitter. Default = 3.
#' @param PSF_diameter The diameter of the PSF in arcseconds. Default = 1".
#' @param SUB Subsampling ratio to be used by the Tiny Tim.
#' @return A bash file with Tiny Tim commands.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitPSFgenerator}}
#' @examples -
#' @export

GRAFitTinyTimACS <- function(wrk_dir = NULL,
                             output_dir = NULL,
                             name,
                             CCDCHIP,
                             x,
                             y,
                             filter = "f814w",
                             focus,
                             ebmv = 0,
                             jitter = 3,
                             PSF_diameter = 1,
                             SUB = 1){

  setwd(output_dir)
  psf_runfile <- paste(output_dir,'/',name,'_run',sep='')
  #     psf_runfile <- paste(name,'_run',sep='')
  fileConn <- file( psf_runfile )
  if (!missing(jitter)) {
    writeLines(paste('~/tinytim-7.5/tiny1 ', name, "_parameter_file ", "ebmv=",ebmv," jitter=",jitter," << ", "EOF",sep=''), fileConn)
  } else{
    writeLines(paste('~/tinytim-7.5/tiny1 ', name, "_parameter_file ", "ebmv=",ebmv," << ", "EOF",sep=''), fileConn)
  }
  close(fileConn)

  cat(15 ,'\n', file = psf_runfile, append = TRUE)         # camera: 15 -> ACS - Wide Field Channel
  cat(CCDCHIP ,'\n', file = psf_runfile, append = TRUE)    # CCDCHIP
  cat(paste(x, y, sep = " ") ,'\n', file = psf_runfile, append = TRUE)        # x-position
  # cat(y ,'\n',file = psf_runfile, append = TRUE)        # y-position
  cat(filter ,'\n', file = psf_runfile, append = TRUE)     # filter
  cat(1 ,'\n', file = psf_runfile, append = TRUE)          # form of object spectrum: 1 -> from list
  cat(13 ,'\n', file = psf_runfile, append = TRUE)         # 13   G8V      1.16   0.75   0.52   0.94   1.40   1.84  ???
  cat(PSF_diameter ,'\n', file = psf_runfile, append = TRUE)     # the PSF diameter; recomended 3.0 arcesecond
  cat(focus ,'\n', file = psf_runfile, append = TRUE)      # the focus (secondary mirror despace)
  cat(paste(name,'_image', sep = '') ,'\n', file = psf_runfile, append = TRUE)
  cat("EOF" ,'\n', file = psf_runfile, append = TRUE)


  system(paste('source ',output_dir,'/',name,'_run',sep=''))
  Sys.sleep(1)
  system(paste('~/tinytim-7.5/tiny2 ', output_dir,'/', name,'_parameter_file', sep = ''))
  Sys.sleep(3)
  system(paste('~/tinytim-7.5/tiny3 ', output_dir,'/', name,'_parameter_file', ' SUB=', SUB, sep = ''))

  system(paste('rm ', output_dir, "/", name,'_image00_psf.fits', sep = ''))
  system(paste('rm ', output_dir, "/", name, '_parameter_file', sep = ''))
  system(paste('rm ', output_dir, "/", name, '_image.tt3', sep = ''))
  # system(paste('rm ', output_dir, "/", name, '_run', sep = ''))

  setwd(wrk_dir)
}

# END


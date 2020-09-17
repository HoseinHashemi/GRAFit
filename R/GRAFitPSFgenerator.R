#' GRAFit: PSF Generator
#'
#' @description This high-level function generates the HST PSF.
#' @param wrk_dir Working directory, where the GRAFit outputs will be saved.
#' @param output_dir The directory to save the PSF.
#' @param HST_Focus_val_file The directory to the file where the HST focus values are stored.
#' @param PSF_name The name of the output PSF.
#' @param target_loc Location; \code{c(RA,DEC)}.
#' @param loc_Unit target_loc unit; \code{"xy"}: CCD coordination or \code{"deg"}: celestial coordination.
#' @param CCDCHIP Which ACS's CCD chip the object is located in. Allowed values: 1 and 2.
#' @param filter What filter is your imaging data. E.g., \code{"f814w"}. If missing it will be read from the \code{header}.
#' @param EXPSTART Time of the start of exposure.
#' @param jitter Telescope's jitter. Default = 3.
#' @param header Header of images.
#' @param ebmv B-V extinction.
#' @param SUB Subsampling ratio to be used by the Tiny Tim.
#' @param PSF_diameter The diameter of the PSF in arcseconds. Default = 1"
#' @param down_sample_factor A factor by which the user wish the Tiny Tim PSF to be resampled.
#' @param verbose Logical; Verbose.
#' @return The image matrix of the PSF.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{TinyTimACSR}}
#' @examples -
#' @export

GRAFitPSFgenerator <- function(wrk_dir = wrk_dir, GRAFitlib = GRAFitlib, output_dir = output_dir, HST_Focus_val_file = NULL,
                               PSF_name = "psf", target_loc, loc_Unit = c("xy","deg"), CCDCHIP, filter = "f814w", EXPSTART, jitter = 3,
                               header = NULL, ebmv = 0, SUB = 1, PSF_diameter = 1, down_sample_factor, verbose = TRUE ) {


  source(paste(GRAFitlib,'/GRAFitTinyTimACS.R',sep=''))

  if (missing(CCDCHIP)) CCDCHIP= as.numeric( header[which(header == "CCDCHIP")+1] )
  if (missing(filter)) filter = as.character( header[which(header == "FILTER2")+1] )
  if (missing(EXPSTART)) EXPSTART = as.numeric( header[which(header == "EXPSTART")+1] )

  # estimate PSF using TinyTim -------------------------------
  # find the location of the object on CCD in pixel ---------
  if (loc_Unit == "deg") {
    xy <- as.integer(magWCSradec2xy(target_loc[1], target_loc[2], header = header))
  } else if (loc_Unit == "xy") {
    xy= target_loc
  }


  # read the focus values from model. Provided by STScI.
  focus_tab <- read.table(file = paste(HST_Focus_val_file,'/Focus.txt', sep = ""), sep = "",
                          col.names = c("Julian date","month","day","year","time","focus_val"))

  focus <- focus_tab[which.min(abs( EXPSTART - focus_tab$Julian.date )), 6]
  focus = focus+0.24   # The correction for ACS/WFC focus. (WFPC2/PC = 0.23, ACS/HRC = -0.25, WFC3/UVIS = -0.24)
  #     The x and y that are converted to pixel by "radec2xy"
  #     should be scaled by the ACS pixel scale which is= 0.05

 # E(B-V) extinction values calculated from Schlegel+98 dust maps and taken from Laigle+2016
  if ( verbose ) cat("Generating PSF ..." ,'\n')

  GRAFitTinyTimACS(wrk_dir = wrk_dir, output_dir = output_dir, name = PSF_name, CCDCHIP = CCDCHIP , xy[1], xy[2], PSF_diameter= PSF_diameter,
                   filter= filter, focus = focus, jitter = jitter, ebmv = ebmv, SUB = SUB )

  psf = readFITS(paste(output_dir,"/",PSF_name,'_image00.fits', sep = ""))$imDat

  if ( !missing(down_sample_factor) ) {
    psf <- profitDownsample(psf, down_sample_factor) }

  return(psf)
}

# END

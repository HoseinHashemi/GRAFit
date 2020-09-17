#' GRAFit: Find the Center of Raw Image Frame.
#'
#' @description This function calculates the center of each image frame and saves coordinates into a file. This is useful for locating each galaxy between large number of frames.
#' @param data_dir The directory where the imaging data (frames) is stored.
#' @return A file containing the coordination of each frame's center.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitFrameCenFinder}}
#' @examples -
#' @export

# This routin reads all the raw frames' headers and then writes out the location (RA & DEC) of the frams' center into a csv file: "frame_cen_cor.csv"

GRAFitFrameCenFinderRawFrame <- function(data_dir = NULL) {

  frame_names_full = list.files(path = data_dir, full.names = TRUE, pattern = "*.fits")
  frame_names = list.files(path = data_dir, full.names = FALSE, pattern = "*.fits")

  RA_ref <- vector()
  DEC_ref <- vector()

  foreach (i=1:length(frame_names)) %do% {
    for( j in seq(2,5, by = 3)  ) {
      # cat(c(i,j), "\n")
      cat('Reading FITS header:: ',frame_names[i], '\n')
      frame_hdr = read.fits(frame_names_full[i])$hdr

      CCDCHIP = as.numeric(frame_hdr[[j]][frame_hdr[[j]][,1]   == "CCDCHIP",2][[1]])
      RA_ref[i] = as.numeric(frame_hdr[[j]][frame_hdr[[j]][,1]   == "CRVAL1",2][[1]])
      DEC_ref[i] = as.numeric(frame_hdr[[j]][frame_hdr[[j]][,1]   == "CRVAL2",2][[1]])

      (trgt_cor <- data.frame(RA_ref = RA_ref[i], DEC_ref = DEC_ref[i], CCDCHIP = CCDCHIP, frame_name = frame_names[i]))
      write.table(trgt_cor, file = paste(data_dir,'/FrameRefRaw.csv',sep = ""), append = T, col.names = F, row.names = F,sep = ",")
    }




  }

  return(file_names)
}

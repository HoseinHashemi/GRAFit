#' GRAFit: Find the Center of Image Frame.
#'
#' @description This function calculates the center of each image frame and saves coordinates into a file. This is useful for locating each galaxy between large number of frames.
#' @param data_dir The directory where the imaging data (frames) is stored.
#' @return A file containing the coordination of each frame's center.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitFrameFinder_v2}}, \code{\link[GRAFit]{GRAFitFrameCenFinderRawFrame}}
#' @examples -
#' @export


# This function reads all the frames' headers and writes out the location (RA & DEC) of the frams' center into a csv file: "frame_cen_cor.csv"

GRAFitFrameCenFinder <- function(data_dir = NULL) {

  file_names = list.files(path = data_dir , full.names = FALSE, pattern = "*.fits")
  RA_TARG <- vector()
  DEC_TARG <- vector()

  foreach (i=1:length(file_names)) %do% {
    cat('Reading FITS header::', file_names[i], '\n')
    hdr <- readFITS(paste(data_dir,'/', file_names[i], sep = ''))$hdr
    RA_TARG[i] = as.numeric(hdr[which(hdr == "CRVAL1") + 1])
    DEC_TARG[i] = as.numeric(hdr[which(hdr == "CRVAL2") + 1])

    (trgt_cor <- data.frame(RA_TARG[i], DEC_TARG[i], file_names[i]))
    write.table(trgt_cor, file = paste(data_dir,'/FrameRef.csv',sep = ""), append = T, col.names = F, row.names = F,sep = ",")

  }

}
# END

#' GRAFit: Find the frame within which the galaxy is located
#'
#' @description This function recieves the location (\code{RA,DEC}) of an object and returns the frame within which the galaxy is located. The frame can then be used for cutouts.
#' @param data_dir The directory where the imaging data (frames) is stored.
#' @param target_loc target location; \code{c(RA,DEC)}.
#' @param FrameType \code{"driz"} for drizzled ACS/HST image or \code{"raw"} for raw ACS/HST image.
#' @return A string specifying the directory/frame.name.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitFrameCenFinder}}, \code{\link[GRAFit]{GRAFitFrameCenFinderRawFrame}}
#' @examples -
#' @export

GRAFitFrameFinder_v2 <- function(data_dir = NULL,
                                 target_loc = NULL,
                                 FrameType = "driz") {

  # ########## READ FRAMES FITS HEADERS ####
  # Read FITS files and write TARGET'S RA AND DEC (i.e. the centre of frames)
  # in a file name: target_cor.csv  ----------------------------

  if (FrameType == "driz") {

    if (!file.exists(paste(data_dir,'/FrameRef.csv',sep = ""))) GRAFitFrameCenFinder( data_dir = data_dir )   # This function reads all the fits header and save them in a file named: target_cor.csv
    file_names = list.files(path = data_dir ,full.names = TRUE, pattern="*.fits") # list all the image (.fits) names in the data_dir.
    FrameRefCor = read.csv(file = paste(data_dir,'/FrameRef.csv',sep = ""),sep = ",")

  } else if (FrameType == "raw") {

    if (!file.exists(paste(data_dir,'/FrameRefRaw.csv',sep = ""))) GRAFitFrameCenFinderRawFrame( data_dir = data_dir )   # This function reads all the fits header and save them in a file named: target_cor.csv
    file_names = list.files(path = data_dir ,full.names = TRUE, pattern="*.fits") # list all the image (.fits) names in the data_dir.
    FrameRefCor = read.csv(file = paste(data_dir,'/FrameRefRaw.csv',sep = ""),sep = ",")
  }

  FrameRefRA <- FrameRefCor[,1]
  FrameRefDEC <- FrameRefCor[,2]
  trgt_cor <- NULL

  ang_sep <- vector()
  # foreach (i=1:length(file_names)) %do% {
  # for( i in 1:length(file_names)) {
  #   # ang_sep[i] = acos( sin(target_loc[2])*sin(FrameRefDEC[i])+cos(target_loc[2])*cos(FrameRefDEC[i])*cos(target_loc[1]-FrameRefRA[i]) )
  #   ang_sep[i] = cenang(a1 = FrameRefRA[i], d1 = FrameRefDEC[i], a2 = target_loc[1], d2 = target_loc[2])
  # }
  #
  # # frame_im_name = file_names[which.min(ang_sep)]
  # frame_im_name = paste(data_dir,"/",FrameRefCor[which.min(ang_sep),3],sep = "")

  frameMatch = coordmatchsing(RAref = target_loc[1], Decref = target_loc[2],
                              coordcompare =  FrameRefCor[,1:2],  rad = 1, radunit = "deg")

  # frame_im_name = paste(data_dir,"/",FrameRefCor[frameMatch$bestmatch[1],3],sep = "")
  frame_im_name = paste(data_dir,"/", FrameRefCor[frameMatch$ID[1], 3], sep = "")
  return(frame_im_name)

}

# END

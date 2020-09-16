# Author: Hosein Hashemizadeh as part of GRAFit package. 
# target_loc: is the Ra and Dec of the target object.
# FrameType: "driz" & "raw" for the drizzled and raw ACS frames, repectively. 

GRAFitFrameFinder <- function(GRAFitlib = GRAFitlib, data_dir = NULL, target_loc = NULL, FrameType = "driz") {

  source(paste(GRAFitlib,'/GRAFitFrameCenFinder.R',sep = ''))
  source(paste(GRAFitlib,'/GRAFitFrameCenFinderRawFrame.R',sep = ''))
  # ########## READ FRAMES FITS HEADERS ####
  # Read FITS files and write the TARGET'S RA AND DEC (i.e. the centre of frames)
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
  for( i in 1:length(file_names)) {
    # ang_sep[i] = acos( sin(target_loc[2])*sin(FrameRefDEC[i])+cos(target_loc[2])*cos(FrameRefDEC[i])*cos(target_loc[1]-FrameRefRA[i]) )
    ang_sep[i] = astro::cenang(a1 = FrameRefRA[i], d1 = FrameRefDEC[i], a2 = target_loc[1], d2 = target_loc[2])
  }
  
  # frame_im_name = file_names[which.min(ang_sep)]
  frame_im_name = paste(data_dir,"/",FrameRefCor[which.min(ang_sep),3],sep = "") 
  return(frame_im_name)
  
}

# END

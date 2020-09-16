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

  # return(file_names)
}


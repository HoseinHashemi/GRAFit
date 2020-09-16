GRAFitPSFgenerator_4stack <- function(wrk_dir = wrk_dir, GRAFitlib = GRAFitlib, output_dir = output_dir, data_dir = data_dir, PSF_name = "psf", 
                               target_loc, loc_Unit = c("xy","deg"), CCDCHIP, filter, EXPSTART, jitter, header=NULL, ebmv=0, tiny_SUB = 1, 
                               down_sample_factor, verbose = TRUE ) {
  
  
  source(paste(GRAFitlib,'/GRAFitTinyTimACS.R',sep=''))
  lckUP_tab1 <- read.csv(file = paste(data_dir,'/look_up_table_EXP1.csv', sep = ""),sep = ",")
  lckUP_tab2 <- read.csv(file = paste(data_dir,'/look_up_table_EXP2.csv', sep = ""),sep = ",")
  lckUP_tab3 <- read.csv(file = paste(data_dir,'/look_up_table_EXP3.csv', sep = ""),sep = ",")
  lckUP_tab4 <- read.csv(file = paste(data_dir,'/look_up_table_EXP4.csv', sep = ""),sep = ",")
  
  CCDCHIP_1 = lckUP_tab1$CCDCHIP[lckUP_tab1$RA==target_loc[1] & lckUP_tab1$DEC==target_loc[2],] 
  
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
  focus_tab <- read.table(file = paste(data_dir,'/Focus.txt', sep = ""), sep = "", 
                          col.names = c("Julian date","month","day","year","time","focus_val"))
  
  focus <- focus_tab[which.min(abs( EXPSTART - focus_tab$Julian.date )), 6]
  focus = focus+0.24   # The correction for ACS/WFC focus. (WFPC2/PC = 0.23, ACS/HRC = -0.25, WFC3/UVIS = -0.24) 
  #     The x and y that are converted to pixel by "radec2xy" 
  #     should be scaled by the ACS pixel scale which is= 0.05
  
  # E(B-V) extinction values calculated from Schlegel+98 dust maps and taken from Laigle+2016
  if ( verbose ) cat("Making PSF ..." ,'\n')
  
  GRAFitTinyTimACS(wrk_dir = wrk_dir, output_dir = output_dir, name = PSF_name, CCDCHIP = CCDCHIP , xy[1], xy[2],
                   filter= filter, focus = focus, jitter = jitter, ebmv = ebmv, SUB = tiny_SUB )
  
  psf = readFITS(paste(output_dir,"/",PSF_name,'_image00.fits', sep = ""))$imDat
  
  if ( !missing(down_sample_factor) ) {
    psf <- profitDownsample(psf, down_sample_factor) }
  
  return(psf)
}

# END

HST_PSF_maker <- function(wrk_dir='~/Desktop/PhD/HST/code/multi_fit', ra=ra,dec=dec) {

  RA_TARG <- trgt_cor[,1]
  DEC_TARG <- trgt_cor[,2]
  trgt_cor<-0
  foreach (j=1:length(ra)) %dopar% {
    
    # input_dir=NULL
    # output_dir=NULL
    ra_deg <- ra[j]
    dec_deg <- dec[j]
    ra_form <- format(ra_deg, digits = 6)
    dec_form <- format(dec_deg, digits = 6)
    system(paste('mkdir ', wrk_dir,'/',ra_form,'+',dec_form, sep=''))
    system(paste('mkdir ', wrk_dir,'/',ra_form,'+',dec_form,'/output', sep=''))
    system(paste('mkdir ', wrk_dir,'/',ra_form,'+',dec_form,'/input', sep=''))
    output_dir <- paste(wrk_dir,'/',ra_form,'+',dec_form,'/output', sep='')
    input_dir  <- paste(wrk_dir,'/',ra_form,'+',dec_form,'/input', sep='')
    
    # estimate PSF using TinyTim -------------------------------
    # find the location of the object on CCD in pixel. ---------s
    xy <- radec2xy(ra_deg, dec_deg, CRVAL1=CRVAL1, CRPIX1 = CRPIX1, 
                   CRVAL2=CRVAL2, CRPIX2=CRPIX2, CD1_1=CD1_1,
                   CD1_2=CD1_2, CD2_1=CD2_1, CD2_2=CD2_2)
    
    # read the focus values from model. Provided by STScI. 
    focus_tab <- read.table(file = paste(data_dir,'/Focus.txt', sep = ""), sep = "", 
                            col.names = c("Julian date","month","day","year","time","focus_val"))
    
    focus <- focus_tab[which.min(abs(EXPSTART-focus_tab$Julian.date)),6]
    
    # The x and y that are converted to pixel by "radec2xy" should be scaled by the ACS pixel scale which is= 0.05
    # TinyTimACSR('psf', CCDCHIP , xy[1]*0.05, xy[2]*0.05, FILTER, focus, PSF_dir = PSF_dir)
    
    
  }     # Loop on all objects. 
  
}
  
  
  
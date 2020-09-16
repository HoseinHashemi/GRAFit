TinyTimACSR <- function(input_dir, name, chip, x, y, filter, focus, exBmV=0, jitter=0){
  setwd(input_dir)
  psf_runfile <- paste(input_dir,'/',name,'_run',sep='')
  #     psf_runfile <- paste(name,'_run',sep='')
  fileConn<-file(psf_runfile)
  writeLines(paste("tiny1 ",name,"_parameter_file ", "ebmv=",exBmV,"jitter=",jitter," << EOF",sep=''), fileConn)
  close(fileConn)
  
  cat(15 ,'\n',file=psf_runfile,append=TRUE)      # camera: 15 -> ACS - Wide Field Channel 
  cat(chip ,'\n',file=psf_runfile,append=TRUE)    # CCDCHIP         
  cat(x ,'\n',file=psf_runfile,append=TRUE)       # x-position      
  cat(y ,'\n',file=psf_runfile,append=TRUE)       # y-position
  cat(filter ,'\n',file=psf_runfile,append=TRUE)  # filter
  cat(1 ,'\n',file=psf_runfile,append=TRUE)       # form of object spectrum: 1 -> from list
  cat(13 ,'\n',file=psf_runfile,append=TRUE)      # 13   G8V      1.16   0.75   0.52   0.94   1.40   1.84  ???
  cat(2.0 ,'\n',file=psf_runfile,append=TRUE)     # the PSF diameter; recomended 3.0 arcesecond
  cat(focus ,'\n',file=psf_runfile,append=TRUE)   # the focus (secondary mirror despace)
  cat(paste(name,'_image',sep='') ,'\n',file=psf_runfile,append=TRUE)
  cat("EOF" ,'\n',file=psf_runfile,append=TRUE)
  
  
  system(paste('source ',input_dir,'/',name,'_run',sep=''))
  Sys.sleep(1)
  system(paste('tiny2 ',input_dir,'/',name,'_parameter_file',sep=''))
  Sys.sleep(3)
  system(paste('tiny3 ',input_dir,'/',name,'_parameter_file',sep=''))
  
  system(paste('rm ', input_dir, '/*_psf.fits',sep= '')) 
  system(paste('rm ', input_dir, '/*_parameter_*',sep= ''))
  system(paste('rm ', input_dir, '/*.tt3',sep= ''))
  # system(paste('rm ', input_dir, '/psf_run',sep= ''))
  # system(paste('mv', 'psf_image*.fits', input_dir,sep= ' '))
  
  setwd(wrk_dir)
}



# GRAFitSBprofile calculates the surface brightness in each individual pixel and make a profile.
# ----------------------
# credit: Hosein Hashemi

GRAFitSBprofile <- function(image=image,main_source=main_src,model=NULL,segim=segim,comp=c("bd","b","d"), centerPos=c(optim_xcen1,optim_ycen1),
                   pixel_scale=0.03, zeropoint = NULL, plot=TRUE, modelOPlot = TRUE, title= NULL,col=NULL,legend=TRUE,
                   legend_lab= c("Data", "Model+Noise", "Bulge", "Disk"),
                   legend_col= c('grey', 'green', 'red', 'blue')) {

  src_pix=image
  src_pix[segim != main_source$segID]=0
  r <- vector()
  mu <- vector()
  rDir <- vector()
  p=0
  
  for(i in 1:nrow(src_pix)) {
    for(j in 1:nrow(src_pix)){
      if (src_pix[i,j] > 0) {
        p=p+1
        r[p]=sqrt((centerPos[1]-i)**2 + (centerPos[2]-j)**2)*pixel_scale
        mu[p]=zeropoint-2.5*log10(src_pix[i,j]/(pixel_scale**2))
        # if(centerPos[1]-i > 0) rDir[p]= -r[p] else rDir[p]= r[p]
      }
    }
  }                                                  
  
  if (modelOPlot) {
    src_pix=model
    segim=segim
    src_pix[segim != main_source$segID]=0
    r_m <- vector()
    mu_m <- vector()
    rDir_m <- vector()
    p=0
    for(i in 1:nrow(src_pix)) {
      for(j in 1:nrow(src_pix)){
        if (src_pix[i,j] > 0) {
          p=p+1
          r_m[p]=sqrt((centerPos[1]-i)**2 + (centerPos[2]-j)**2)*pixel_scale    
          mu_m[p]=zeropoint-2.5*log10(src_pix[i,j]/(pixel_scale**2))
          # if(centerPos[1]-i > 0) rDir_m[p]= -r[p] else rDir_m[p]= r[p]
        }
      }
    }
    #par(mar = rep(5, 4))
    
    if(plot) {
      if(comp=="bd") {
        magplot(r, mu, xlim = c(min(r),max(r)), ylim = c(27,min(mu)), 
                pch=19, cex=0.3, xlab = "Projected Distance to Center /asec" , 
                ylab = expression(mu(mag/asec^2)), cex.lab=1.5, labels =1, 
                grid = TRUE,main=title,col='grey')        
      }
      
      par(new=TRUE)
      magplot(r_m, mu_m, xlim = c(min(r),max(r)), ylim = c(27,min(mu)), 
              pch=19, cex=0.3, col=col, xlab = "", ylab = "", labels =F)
      if (legend) {
        if (comp=="bd") {
          legend("topright", c(legend_lab[1],legend_lab[2],legend_lab[3],legend_lab[4]),
                 col = c(legend_col[1],legend_col[2],legend_col[3],legend_col[4]), 
                 pch = 19, y.intersp = 1, bg = 'white') 
        } else if (comp=="b" | comp=="d") {
          legend("topright", c(legend_lab[1],legend_lab[2]),col = c(legend_col[1],legend_col[2]), 
                 pch = 19, y.intersp = 1, bg = 'white') 
        }
        
      }
      # return(list(rSrc=rDir, muSrc=mu,rModel=rDir_m,muModel=mu_m))
      # return(list(pl))
    } #else return(list(rSrc=rDir, muSrc=mu,rModel=rDir_m,muModel=mu_m))
    
  } else {
    if(plot) {
      par(mar = rep(5, 4))
      magplot(r, mu, ylim = c(max(mu),min(mu)), pch=19, cex=0.3, xlab = "Projected Distance to Center /asec" ,
              ylab = expression(mu(mag/asec^2)), cex.lab=1.5 ,side = c(1,2,3,4),main=title)  
      legend("topright", c("Data"), pch = 19, bg = 'white')
    }
    # return(list(rSrc=rDir, muSrc=mu))
    
  }

}
# END


GRAFitSBdistr <- function(image=image,main_source=main_src,model=NULL,bulge_Model=NULL,disk_Model=NULL,segim=segim,comp=c("bd","b","d")) {

  src_mdlSB=GRAFitSBprofile(image=image,main_source=main_src,segim=segim,model=LDmodel_noise,comp="bd",
                            col = 'green', modelOPlot=TRUE, legend=TRUE, main = "Surface Brightness Profile")
  
  if (nComp == 2 ) {
    par(new=TRUE)
    B_mdlSB=GRAFitSBprofile(image=image,main_source=main_src,segim=segim,model=bulge_model$z+noise,comp = "b",
                            modelOPlot=TRUE,col='red',legend = FALSE, plot = FALSE)
    D_mdlSB=GRAFitSBprofile(image=image,main_source=main_src,segim=segim,model=disk_model$z+noise,comp = "d",
                            modelOPlot=TRUE,col='blue',legend = FALSE, plot = FALSE)
  }
  spar=0.7;bins=300
  plot.new()
  mean_src <- magrun(src_mdlSB$rSrc*PIXELSCALE,src_mdlSB$muSrc,type = 'mean',bins = bins)
  magplot(mean_src$x,mean_src$yquan[,1], ylim = c(26,18),xlim = c(min(src_mdlSB$rSrc*PIXELSCALE),max(src_mdlSB$rSrc*PIXELSCALE)),col='grey',pch=16, xlab = "Distance to center / asec", ylab = expression(mu(mag/asec^2) ))
  smoothingSpline = smooth.spline(mean_src$x,mean_src$yquan[,1], spar=spar)
  lines(smoothingSpline)
  
  par(new=TRUE)
  mean_model <- magrun(src_mdlSB$rModel*PIXELSCALE,src_mdlSB$muModel,type = 'mean',bins = bins)
  # magplot(mean_model$x,mean_model$yquan[,1], ylim = c(26,18),xlim = c(-3,3),col='green')
  # lines(mean_model$x,mean_model$yquan[,1], ylim = c(26,18),xlim = c(-3,3),col='green')
  smoothingSpline = smooth.spline(mean_model$x,mean_model$yquan[,1], spar=spar)
  lines(smoothingSpline, ylim = c(26,18),xlim = c(-3,3),col='green')
  
  par(new=TRUE)
  mean_B <- magrun(B_mdlSB$rModel*PIXELSCALE,B_mdlSB$muModel,type = 'mean',bins = bins)
  # magplot(mean_B$x,mean_B$yquan[,1], ylim = c(26,18),xlim = c(-3,3),col='red')
  # lines(mean_B$x,mean_B$yquan[,1], ylim = c(26,18),xlim = c(-3,3),col='red')
  smoothingSpline = smooth.spline(mean_B$x,mean_B$yquan[,1], spar=spar)
  lines(smoothingSpline, ylim = c(26,18),xlim = c(-3,3),col='red')
  
  par(new=TRUE)
  mean_D <- magrun(D_mdlSB$rModel*PIXELSCALE,D_mdlSB$muModel,type = 'mean',bins = bins)
  # magplot(mean_D$x,mean_D$yquan[,1], ylim = c(26,18),xlim = c(-3,3),col='blue')
  # lines(mean_D$x,mean_D$yquan[,1], ylim = c(26,18),xlim = c(-3,3),col='blue')
  smoothingSpline = smooth.spline(mean_D$x,mean_D$yquan[,1], spar=spar)
  lines(smoothingSpline, ylim = c(26,18),xlim = c(-3,3),col='blue')
  
  legend("topright", c("Object","Model","Disk","Bulge"),col = c('black','green','blue','red'), 
         lty=1, y.intersp = 1, bg = 'white')
  
  
  
}


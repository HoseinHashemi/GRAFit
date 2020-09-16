# GRAFitAperturePlot function overplots the extracted source apertures on the image. 

GRAFitAperturePlot <- function(x=x, y=y, a=a, b=b, angle=angle, border= 'red') {
  points(x,y,pch= 4, col=border)
  draw.ellipse(x,y,a,b,angle, border = border)
}
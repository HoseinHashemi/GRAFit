# Author: Aaron Robotha as part of magicaxis package.
# Modified by : Hosein Hashemi as part of GRAFit package.

GRAFitri <- function (chains, samples, samptype = "end", grid = FALSE, tick = FALSE) 
{
  
  chains = as.data.frame(chains)
  chaincolnames = colnames(chains)
  Nsamp = dim(chains)[1]
  Npar = dim(chains)[2]
  if (Npar <= 1) {
    stop("Need 2+ parameters!")
  }
  if (missing(samples)) {
    samples = Nsamp
  }
  if (samples > Nsamp) {
    samples = Nsamp
  }
  layout(matrix(1:Npar^2, Npar, Npar)[Npar:1, ])
  meanvec = {
  }
  sdvec = {
  }
  if (samptype == "end") {
    usesamps = (Nsamp - samples + 1):Nsamp
  }
  if (samptype == "ran") {
    usesamps = sample(Nsamp, samples)
  }
  for (i in 1:Npar) {
    meanvec = c(meanvec, mean(chains[usesamps, i]))
    sdvec = c(sdvec, sd(chains[usesamps, i]))
  }
  par(oma = c(8, 8, 3, 3), cex.axis = 1.)
  for (i in 1:Npar) {
    for (j in 1:Npar) {
      par(mar = c(0, 0, 0, 0) )
      xrange = range(chains[usesamps, i])
      yrange = range(chains[usesamps, j])
      if (xrange[1] == xrange[2]) {
        val = xrange[1]
        xrange[1] = val - 0.05
        xrange[2] = val + 0.05
      }
      if (yrange[1] == yrange[2]) {
        val = yrange[1]
        yrange[1] = val - 0.05
        yrange[2] = val + 0.05
      }
      if (i == j) {
        xtemp = chains[usesamps, i]
        if (sd(xtemp) == 0) {
          xtemp = xtemp + rnorm(samples, sd = 0.001)
        }

        plot(density(xtemp), axes = FALSE, main = "", 
             xlim = xrange, lwd = 4, col = c("purple"))
        magaxis(1, grid = grid, grid.col = "lightgrey", 
                labels = FALSE, tick = tick, cex.lab = 2.2)
        abline(v = meanvec[i], lty = 2, col = "red", lwd = 2)
        abline(v = meanvec[i] - sdvec[i], lty = 3, col = "red", lwd = 2)
        abline(v = meanvec[i] + sdvec[i], lty = 3, col = "red", lwd = 2)
        box()
        if (i == 1) {                                                          
          plot.window(xlim = xrange, ylim = yrange)
          magaxis(1:2, xlab = chaincolnames[i], ylab = chaincolnames[j], cex.lab = 2.5, mtline = 4)
        }
      }
      else {
        if (i > j) {
          plot.new()
          plot.window(xlim = xrange, ylim = yrange)
          xtemp = chains[usesamps, i]
          ytemp = chains[usesamps, j]
          if (sd(xtemp) == 0) {
            xtemp = xtemp + rnorm(samples, sd = 0.001)
          }
          if (sd(ytemp) == 0) {
            ytemp = ytemp + rnorm(samples, sd = 0.001)
          }
          magaxis(1:2, grid = grid, grid.col = "lightgrey", 
                  labels = FALSE, tick = tick)
          magcon(xtemp, ytemp, dobar = FALSE, doim = FALSE, 
                 add = TRUE, lty = c(2, 1, 3), xlim = xrange, 
                 ylim = yrange, lwd = 4, col = c("turquoise2"))
          points(meanvec[i], meanvec[j], col = "red", 
                 pch = 4, cex = 6, lwd = 2)
          box()
          abline(v = meanvec[i], lty = 2, col = "red", lwd = 2)
          abline(v = meanvec[i] - sdvec[i], lty = 3, 
                 col = "red", lwd = 2)
          abline(v = meanvec[i] + sdvec[i], lty = 3, 
                 col = "red", lwd = 2)
          if (j == 1) {
            magaxis(1, xlab = chaincolnames[i], cex.lab = 2.5, mtline = 4)
          }
        }
        else {
          plot.new()
          
          plot.window(xlim = xrange, ylim = yrange)
          magaxis(1:2, grid = grid, grid.col = "lightgrey", 
                  labels = FALSE, tick = tick)
          points(chains[usesamps, c(i, j)], pch = 16, 
                  cex = .7, col = c("gold"))
          points(meanvec[i], meanvec[j], col = "red", 
                 pch = 4, cex = 6, lwd = 2)
          box()
          if (i == 1) {
            magaxis(2, ylab = chaincolnames[j], cex.lab = 2.5, mtline = 4)
          }
        }
      }
    }
  }
  output = cbind(mean = meanvec, sd = sdvec)
  rownames(output) = chaincolnames
  return = output
}


#------------------------------

GRAFitPriorsPlot <- function (Data = Data, sigmas = c(2, 2, 2, 2, 5, 5, 0.3, 0.3, 0.3, 0.1, 
                                              30, 30, 0.3, 0.3, Inf, Inf), allowflat = TRUE) 
{
  
  parm = Data$init
  remakeout = profitRemakeModellist(parm = parm, Data = Data)
  new = remakeout$modellist
  parm = remakeout$parm
  init = Data$modellist
  means = unlist(init)
  
  tolog = unlist(Data$tolog)
  tofit = unlist(Data$tofit)
  intervals <- unlist(Data$intervals)
  par(mfrow = c(3,4))
  LL = 0
  parms = unlist(new)
  if (!is.null(tofit)) {
    ps = which(tofit)
  }   else ps = 1:length(parms)

  for (p in 1:14) {
    means = unlist(init)
    xlab = names(unlist(init))

    # main = c("xcen1", "xcen2", "ycen1", "ycen2", "mag1", "mag2",
    #          "re1", "re2", "nser1", "nser2", "ang2", "axrat2")
    
    print(p)
    xlab = xlab[p]
    sd = sigmas[p]
    mean = means[[p]]
    main = substitute( paste( sigma, " = " , sd), list(sd = sd) )
    
    if (!(allowflat && (sigmas[p] == Inf))) {
      parm = parms[[p]]
      mean = means[[p]]
      
      if (tolog[p]) {
        parm = log10(parm)
        mean = log10(mean)
        xlab = paste("log10(",xlab, ")", sep = '')
        main = substitute( paste( sigma, " = " , sd, " dex"), list(sd = sd) )
      }
      
      x <- seq( mean-(2*sd), mean+(2*sd), .01)

      pr = dnorm(x, mean, sigmas[p], log = TRUE)
      if (tofit[p]){
        if ( p < 5 ) {
          magplot(x, pr, log = '', type = 'l', lwd = 2, xlab = xlab, 
                  main = main, cex.main = 2,
                  xlim = c(mean-.05*mean, mean+.05*mean), cex.lab = 1.5)   
          abline(v = c(mean,mean-sd,mean+sd), col = 'red', lty = c(1,2,2), lwd = c(2,1,1))
          
        } else {
          magplot(x, pr, log = '', type = 'l', lwd = 2, xlab = xlab, 
                  main = main, cex.main = 2,
                  xlim = c(mean-2*mean, mean+2*mean), cex.lab = 1.5)   
          abline(v = c(mean,mean-sd,mean+sd), col = 'red', lty = c(1,2,2), lwd = c(2,1,1))
          
        }
      }
       
      # if (tolog[p]) {
      #   magplot(x, pr, log = 'y', type = 'l', lwd = 1.5, main = main) 
      # } else {
        
      # }
      
      
    } else if (allowflat && (sigmas[p] == Inf)) {
      # magplot(x, pr, log = '', type = n, main = main)
      # abline(h = mean, col = 'red', lty = 2, lwd = 1.5)
    }
    
  }

}



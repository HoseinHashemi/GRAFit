# This tool calculates the central surface brightness from Graham & Driver 2005.
# Author: Hosein Hashemizadeh

GRAFitMu0 <- function(n, mag, Re, PixScale) {
  
  if (!missing(PixScale)) {
    Re = Re * PixScale
  } else Re = Re
  bn = 1.9992 * n - 0.3271

  mu0 = mag + 5 * log10(Re) - 2.5 * log10( bn ** (2*n) / (pi * gamma(2 * n + 1))  )
  
  return(mu0)
}

#END
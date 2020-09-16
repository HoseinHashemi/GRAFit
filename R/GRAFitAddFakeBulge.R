GRAFitAddFakeBulge <- function(model = NULL, zeropoint = 34.0) # Add a zero-point magnitude bulge
{
  # <param: model [list]> - The modellist from profitMakeModel()
  # <return: pseudo [list]> - A (pseudo) modellist containing a second component; i.e., a duplicate with "zero magnitude".

  pseudo = model # copy model list
  for (key in names(model$sersic)){pseudo$sersic[[key]][2] = model$sersic[[key]][1]} # duplicate a second component
  pseudo$sersic$mag[1] = zeropoint # set magnitude of bulge (loc=1) to ZERO_POINT
  return(pseudo)
}

# END

# Authot: Hosein Hashemi June 2018.
# This function classify galaxies in different classes regarding their 1D light profile.

GRAFitfitClass <- function(SB = NULL, SBlim = 26) {

  if (is.null(SB)) {
    fitClass = NA
    warning("WARNING: Surface brightness info is not provided, fitClass=NA ")
  } else {

    ellipseSBlim = SB[ SB$total_SB < SBlim, ]
    # DB_ratio <- as.numeric(format(ellipseSBlim$disk_SB/ellipseSBlim$bulge_SB, digits = 5))
    DB_ratio <- ellipseSBlim$disk_SB/ellipseSBlim$bulge_SB
    DB_diff <- ellipseSBlim$disk_SB-ellipseSBlim$bulge_SB
    # interc <- which.min(abs(DB_ratio-1))
    # # interc <- DB_ratio[DB_ratio == 1.000]
    # intercNo = length(interc)
    
    if ( first(DB_diff)/last(DB_diff) < 0 & 
        first(ellipseSBlim$bulge_SB) < first(ellipseSBlim$disk_SB) ) {
      
      fitClass = 0      # classic fit
      
    } else if (first(DB_diff)/last(DB_diff) < 0  & 
               last(ellipseSBlim$bulge_SB) < last(ellipseSBlim$disk_SB)) {
      
      fitClass = 1     # bulge fitting the outer halo or sky
      
    } else if (first(DB_diff)/last(DB_diff) < 0 &
               first(ellipseSBlim$bulge_SB) > first(ellipseSBlim$disk_SB)) {
      
      fitClass = 2     # bulge fitting truncated disk/outer region
      
    } else if (first(DB_diff)/last(DB_diff) > 0 & 
               !all(DB_ratio > 1.) & 
               first(ellipseSBlim$bulge_SB) < first(ellipseSBlim$disk_SB) &
               last(ellipseSBlim$bulge_SB) < last(ellipseSBlim$disk_SB) ){
      
      fitClass = 3     # single component ellipticals or bulge trying to fitt outer halo (B/T > 0.5, nb > 1.5) 
      
    } else if ( first(DB_diff)/last(DB_diff) > 0 & 
                !all(DB_ratio > 1.) & !all(DB_ratio < 1.) &
                first(ellipseSBlim$bulge_SB) > first(ellipseSBlim$disk_SB) &
                last(ellipseSBlim$bulge_SB) > last(ellipseSBlim$disk_SB) ){
      
      fitClass = 4     # disk with perturbations or feature of irregular galaxies 
      
    } else if (first(DB_diff)/last(DB_diff) > 0 &
               all(DB_ratio < 1.) ) {
      
      fitClass = 5  # disk dominated (B/T < 0.5) or pure disk (B/T = 0). Likely single sersic.
      
    } else if (first(DB_diff)/last(DB_diff) > 0 & 
               all(DB_ratio > 1.) ) {
      
      fitClass = 6  # bulge dominated (B/T > 0.5) or pure bulge (B/T = 1). Likely elliptical single sersic.             
      
    } else fitClass = NA 
    
  }
  
return(fitClass = fitClass)
    
}
# DONE!


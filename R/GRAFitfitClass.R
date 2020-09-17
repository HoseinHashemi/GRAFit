#' GRAFit: Fit Class.
#'
#' @description This function returns an integer indicating the fit class according to how many times the 1D profiles of the data and model cross.
#' @param SB A vector of surface brightness values.
#' @param SBlim 5 sigma surface brightness limit of the data typically in units of mag/asec^2. Default = 26.
#' @return An integer between 0-6. 0: normal fit, 1: bulge is likely fitting the outer halo or sky, 2: bulge is likely fitting truncated disk/outer region, 3: single component ellipticals or bulge trying to fit outer halo, 4: likely a clumpy disk with perturbations or strong feature of irregular galaxies, 5: likely a disk dominated (B/T < 0.5) or pure disk (B/T = 0) system. 6: likely a bulge dominated or pure bulge (B/T = 1) system. Likely an elliptical single sersic.
#' @author Hosein Hashemizadeh
#' @seealso \code{\link[GRAFit]{GRAFitEllipsePlot}}
#' @examples -
#' @export

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

      fitClass = 3     # single component ellipticals or bulge trying to fit outer halo (B/T > 0.5, nb > 1.5)

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
# END!


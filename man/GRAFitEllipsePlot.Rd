% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRAFitEllipsePlot.R
\name{GRAFitEllipsePlot}
\alias{GRAFitEllipsePlot}
\title{GRAFit: Plot Isophotal Surface Brightness for Pseudo-Ellipses}
\usage{
GRAFitEllipsePlot(
  Data,
  modellist,
  bulgeloc = 1,
  diskloc = 2,
  pixscale = 1,
  FWHM = 0.09,
  SBlim = 26,
  df = 100,
  raw = FALSE,
  plot = TRUE,
  main_src = NULL
)
}
\arguments{
\item{Data}{Data of class profit.data. The standard Data structure containing all inputs, images, modelsetup etc.}

\item{modellist}{Model list (see \code{\link[ProFit]{profitMakeModel}} for more details).}

\item{bulgeloc}{Location ID of bulge component in the Sersic list provided in modellist}

\item{diskloc}{Location ID of disk component in the Sersic list provided in modellist.}

\item{FWHM}{The full width half max of the PSF in units of arc seconds. A vertical line is drawn at half this number (since we are plotting radius). The fits inside of the region is inherently hard because it is well within the PSF convolution kernel}

\item{SBlim}{5 sigma surface brightness limit of the data typically in units of mag/asec^2. Default = 26.}

\item{df}{Degrees of freedom to use for spline fitting. Lower if the lines look too wavy.}

\item{raw}{Logical; if FALSE (the default) then a smooth spline is used to represent the data and model 1D profiles. This smooths out deprojection noise caused by the PSF often being non-smooth. If TRUE then the raw pixel surface brightness values are shown. These will show much more scatter, but the trends ought to be very similar. If the raw and smooth 1D plots differ significantly then the df flag probably needs to be changed to improve the smoothing. Notice that when the raw pixel values are plotted the shaded error polygon is very hard to see (it is usually subdominant compared to the pixel scatter created during deprojection that has both the Normal pixel error and the PSF induced deprojection error).}

\item{plot}{Logical; should the plot be generated.}

\item{main_src}{A list containing the specification of the main source in the cutout. This should be generated by the \code{ProFound}.}
}
\value{
The 1D flux density profile plus a list containing the surface brightness of each pixel and fit class, i.e., how the bulge and disk profiles cross. See \code{\link[GRAFit]{GRAFitfitClass}} for more details.
}
\description{
This function is a modification of \code{profitEllipsePlot} to plot the 1D profile of the flux intensity.
}
\examples{
-
}
\seealso{
\code{\link[GRAFit]{GRAFitfitClass}}, \code{\link[ProFit]{profitEllipsePlot}}
}
\author{
Hosein Hashemizadeh
}

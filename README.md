
<!-- README.md is generated from README.Rmd. Please edit that file -->
GRAFit (R package)
==================

<!-- badges: start -->
<!-- badges: end -->
GRAFit is an modular automatic pipeline for galaxy structural decomposition originally designed for the Hubble Space Telescope (HST) imaging data, the Advanced Camera for Surveys (ACS) in particular. However, it can be used for any standard astronomical imaging data. <!-- At it's core GRAFit uses ProFit galaxy fitting package. -->

Installation
------------

### Install R

Get R for OSX (Mac), Linux and Windows from: <https://cloud.r-project.org>

### Install ProFit

GRAFit uses ProFit, a Bayesian galaxy profile fitting tool.

Get ProFit with:

``` r
install.packages('devtools')
library(devtools)
install_github("ICRAR/ProFit")
library(ProFit)
```

See [ProFit](https://github.com/ICRAR/ProFit) if you faced any issues with installing ProFit or it's dependencies.

### Install ProFound

GRAFit uses ProFound for photometry and source extraction.

Get ProFound with:

``` r
install.packages('devtools')
devtools::install_github("asgr/ProFound")
library(ProFound)
```

This should work properly, and the dependencies should be installed automatically.

See [ProFound](https://github.com/asgr/ProFound) if you faced any issues with installing ProFound.

### Install Tiny Tim

GRAFit uses Tiny Tim for generating HST Point Spread Function (PSF).

Get Tiny Tim from [here](https://www.stsci.edu/hst/instrumentation/focus-and-pointing/focus/tiny-tim-hst-psf-modeling)

### Install GRAFit

<!-- You can install the released version of GRAFit from [CRAN](https://CRAN.R-project.org) with:-->
<!-- ``` r -->
<!-- install.packages("GRAFit") -->
<!-- ``` -->
Install GRAFit from [GitHub](https://github.com/hoseinhashemi/GRAFit) with:

``` r
<!-- devtools::install_github("hoseinhashemi/GRAFit") -->
# library(GRAFit)
```

Example
-------

In principle, you only need to use one function of GRAFit; `GRAFitMaster`. This function calls whatever it needs to do the fitting.

You should provide: - a working directory where the outputs will be saved.

-   a directory where the your imaging data is stored.

-   If any, a directory where you have stored your pre-generated Point Spread Function (PSF) for each object. Otherwise, you may leave a PSF in your wrk\_dir or alternatively GRAFit will generate a PSF using Tiny Tim.

This is a basic example which shows you how to run GRAFit:

``` r
# GRAFitMaster(wrk_dir = "wrk_dir", 
#              data_dir = "data_dir", 
#              PSF_dir = "PSF_dir", 
#              threadMode = 0, 
#              ncores = 1, 
#              nComp= 2, 
#              optimMode = 'MCMC', 
#              object_list = object_list)
```

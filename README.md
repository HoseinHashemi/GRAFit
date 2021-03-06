
<!-- README.md is generated from README.Rmd. Please edit that file -->
GRAFit (R package)
==================

<!-- badges: start -->
<!-- badges: end -->
**GRAFit** is an automatic modular pipeline for galaxy structural decomposition (also known as bulge-disk decomposition) originally designed for the Hubble Space Telescope (HST) imaging data, the Advanced Camera for Surveys (ACS) in particular. However, it can be used for any standard astronomical imaging data. **GRAFit** is able to fit galaxies with different combination of models including double Sersic profile and single Sersic profile. <!-- At it's core GRAFit uses ProFit galaxy fitting package. -->

Installation
------------

### Install R

Get **R** for OSX (Mac), Linux and Windows from: <https://cloud.r-project.org>

### Install ProFit

For fitting process **GRAFit** uses **ProFit**, a Bayesian galaxy profile fitting tool.

Get **ProFit** with:

``` r
install.packages('devtools')
library(devtools)
install_github("ICRAR/ProFit")
library(ProFit)
```

See [ProFit](https://github.com/ICRAR/ProFit) if you faced any issues with installation and/or dependencies.

### Install ProFound

Before we model a galaxy we need an approximate estimation of some photometrical parameters such as magnitude, central positions, position angle, etc. **GRAFit** will parse these values to **ProFit** as initial inputs of our model. **GRAFit** uses **ProFound** for photometry, source extraction and generating segmentation and mask maps.

Get **ProFound** with:

``` r
devtools::install_github("asgr/ProFound")
library(ProFound)
```

This should work properly, and the dependencies should be installed automatically.

See [ProFound](https://github.com/asgr/ProFound) if you faced any issues with installation.

### Install Tiny Tim

To simulate the optical effects on our model we need to convolve the with a Point Spread Function (PSF). **GRAFit** uses **Tiny Tim** for generating HST PSF.

Get **Tiny Tim** from [here](https://www.stsci.edu/hst/instrumentation/focus-and-pointing/focus/tiny-tim-hst-psf-modeling)

The user can also provide the PSF in their working directory or any other directories and provide the direction to them. This way **GRAFit** will prioretize the provided PSFs and will not generate new PSFs anymore.

### Install GRAFit

<!-- You can install the released version of GRAFit from [CRAN](https://CRAN.R-project.org) with:-->
<!-- ``` r -->
<!-- install.packages("GRAFit") -->
<!-- ``` -->
You are now ready to install **GRAFit** from [GitHub](https://github.com/hoseinhashemi/GRAFit) with:

``` r
devtools::install_github("hoseinhashemi/GRAFit")
library(GRAFit)
```

### Dependencies

Most of the required dependencies should be installed by this stage. In case any of the dependencies are failed to install, you can install them with.

``` r
install.packages(c('fftw', 'R2Cuba', 'RColorBrewer', 'EBImage', 
                   'LaplacesDemon', 'imager', 'magicaxis', 
                   'FITSio', 'data.table'))
```

The segmentation map required for **GRAFit** is generated by **ProFound** using **EBImage** package. If this package was not installed by this stage install it with:

``` r
install.packages("BiocManager")
BiocManager::install("EBImage")
```

### Modules

**GRAFit** comes with a number of modules. For more details about each module and get the documentation for functions run the help operator in R environment with `?functionname` or `help(functionname)`. For example:

``` r
?GRAFitDynamo_v2
help(GRAFitDynamo_v2)
```

Example
-------

**GRAFit** is originally designed for a large number of galaxies distributed between one or more imaging frames, although there is no restriction for fitting one galaxy in a single cutout. In principle, you need to interact with only one **GRAFit** function; `GRAFitMaster`. This function calls all other modules and dependencies when/where it needs during the process.

You should provide:

-   a working directory where the outputs will be saved.

-   a directory where your imaging data is stored.

-   If any, a directory where you have stored your pre-generated Point Spread Function (PSF) for each object. Otherwise, you may leave a PSF in your `wrk_dir` or alternatively GRAFit will generate a HST PSF using **Tiny Tim**. In the later you will beed to provide the PSF specifications required by **Tiny Tim**.

For a simple run, I provide a folder including a COSMOS HST/ACS frame together with the associated PSF and an input catalogue consisting of a random galaxy. Download the whole folder from [here](https://www.dropbox.com/sh/pb8xri702tz03ka/AAB7hwuy89jhyqRHUwuMrBK0a?dl=0) and put it where you wish in you local directory, suppose in your Desktop.

Then run the code below for a basic example of **GRAFit** fitting.

``` r

inCat = read.csv('~/Desktop/GRAFit_test_Run/Input.csv')

GRAFitMaster(wrk_dir = "~/Desktop/GRAFit_test_Run/",
             data_dir = "~/Desktop/GRAFit_test_Run/",
             threadMode = 0,
             ncores = 1,
             nComp= 2,
             optimMode = 'optim',
             object_list = inCat)
```

-   Note that `object_list` should include at least the RA and DEC of your target galaxy.

This run uses a relatively simple `optim` optimizer to save you some time. More robust optimization that **GRAFit** provides is MCMC for which you need to specify optimMode = 'MCMC'. See `GRAFitMaster` documentation for more details and options.

This run will genrate a folder named with the galaxy ID (or RA/DEC if ID is not provided in inCat) containing all outputs associated with this galaxy such as figures of initial and optimized models, 1D profiles etc. In addition it will generate a master catalogue (MasterCat.csv) in you `wrk_dir` directory including all fitting parameters as well as a ProCat.csv including the results of photometric measurements by **ProFound**. Another directory could be also generated containing work spaces for which you need to specify `keep_wrk_space = TRUE`.

Contributor(s)
--------------

Hosein Hashemizadeh

License
-------

L-GPL3

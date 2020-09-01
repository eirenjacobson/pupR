<img src="vignettes/figures/Logos.png" align="right" height="50px"/>

# The pupR package

#### T. A. Øigård and M. Biuw

# Overview
The `pupR` package is a R-package for estimating the pup production of harp and hooded seals in the West Ice (along the coast of Greenland). Abundance estimation of other populations is possible as long as similar input data is available. The package fits available (photographic) count data (from aerial photos) and estimates the abundance of the population of interest. Any type of count data can be used to estimate the abundance or total numbers for a specified area.

Instructions on how to use `pupR` is found in the vignette.


# Installation

Dependencies of other R-packages will be automatically installed when installing the `pupR` package.

The most recent version of `pupR` is hosted on a git repository at
<https://github.com/NorskRegnesentral/pupR.git>.


To install the R-package directly from the repository use the following command (note: the R-package devtools has to be installed first)
``` r
devtools::install_github("https://github.com/NorskRegnesentral/pupR.git", build_vignettes = TRUE)
``` 

In order to load the package type:
```{r}
library(pupR)
```

Instructions on how to use `pupR` is found in the vignette.

To load the vignette type:
```{r}
vignette("howToUse_pupR",package = "pupR")
```

The vignette will be opened in a system specified viewer. For a HTML version of the vignette (recommended) use

```{r}
browseVignettes()
```
and scroll down to the `pupR` package.

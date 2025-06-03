First it's necessary to have devtools installed:

```r
# Installing the devtools
install.packages("devtools")

library(devtools)
```
It's required to have Rtools to install the package. It can be downloaded at: https://cran.r-project.org/bin/windows/Rtools/

Then the user must install de package using devtools and the appropriate path:

```r
devtools::install_github("ThiagoMCarneiro/GeoPoisson")
library(GeoPoisson)
```

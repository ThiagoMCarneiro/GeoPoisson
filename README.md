This package is based on the non-Homogeneous Poisson model by Morales et al. (2016) that can be found at: https://repositorio.ufrn.br/items/70c7e1a7-ea8a-41bc-9034-93920cbb151e

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
Use the packages below:

```r
library(MASS)
library(Matrix)
```

Finally, to use the main function:

```r
result <- GeoPoisson(data,limuser,0.001,0.001,n_iter,burn_in,X,loca)

```
The variable "data" are the observations from stations, limuser is a threshold to count extreme data, c1 and d1 are hyperparameters (0.001 is recommended for both), n_iter is the number of iterations, burn_in is the number of iterations that will be ignores, X and loca are the matrices of the location of stations (and X has a 1 coefficient for the intercept).

The variable "result" is a list of 6 variables, result[[1]] to result [[6]] and it will have (n_iter - burn_in) number of observations of alpha, beta, b, v, W (a matrix) and Psi (a matrix).

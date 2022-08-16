
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MPEwSVMs

<!-- badges: start -->
<!-- badges: end -->

The goal of MPEwSVMs is to perform the multiclass probability estimation
with weighted SVMs, using pairwise coupling, OVA and linear time
algorithm from: “Zeng, L. and Zhang, H. H. (2022). Linear algorithms for
nonparametric multiclass probability estimation”. \[link\]
(<https://arxiv.org/abs/2205.12460>).

## Installation

You can install the development version of MPEwSVMs from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("zly555/MPEwSVMs")
```

## Example

Here we use a example to demonstrate the multiclass probability
estimation with “MPEwSVMs” package. We use a three class non-linear
simulation example from [link](https://arxiv.org/abs/2205.12460v1)of
example 5 (page 18).

### Install the required library

``` r
library(quadprog)
library(data.table, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(tictoc, warn.conflicts = FALSE)
library(MPEwSVMs)
```

### The function to generate the simulation data

``` r
generator_ex2 <- function(n, seed=1){
  
  set.seed (seed)

  x1 <- runif(n, min=-3, max=3)
  x2 <- runif(n, min=-6, max=6)

  x <- cbind(x1,x2)

  data.dim <- dim(x)[2]

  f1 <- function(x) {return(-x[1] +0.1*x[1]^2-0.05*x[2]^2+0.1)}
  f2 <- function(x) {return(-0.2*x[1]^2 +0.1*x[2]^2-0.2)}
  f3 <- function(x) {return(x[1]+0.1*x[1]^2-0.05*x[2]^2+0.1)}

  # calculate the pjx for each class
  pyx <- matrix(0, nrow=n, ncol =3)

  for (i in 1:3){
    pyx[,i] <- apply(x, 1, function(x) {exp(do.call(paste0("f",i), list(x)))/(exp(f1(x))+exp(f2(x))+exp(f3(x)))})
  }

  # generate moltinorminal random variable for y label, 1 indicate pick that class label
  y <- apply(pyx, 1, function(x) rmultinom(n=1, size =1, prob=x))

  # generate y label {1,2,3}
  y <- apply(y*c(1,2,3), 2, function(x) x[x > 0])

  data <- data.frame(cbind(x, y, pyx))
  colnames(data) <- c('x1','x2','ylabel','p1x','p2x','p3x')

  return(list(data = data, data_dim = data.dim))

}
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

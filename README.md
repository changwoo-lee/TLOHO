
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TLOHO

<img src="TLoHo_longvertical_text.png" width="200">
<!-- badges: start --> <!-- badges: end -->

TLOHO is an accompanying R package of the paper:

Lee, C. J., Luo, Z. T., & Sang, H. (2021). T-LoHo: A Bayesian
Regularization Model for Structured Sparsity and Smoothness on Graphs.
*Advances in Neural Information Processing Systems 34 (NeurIPS 2021)*

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("changwoo-lee/TLOHO")
```

## Example

### Running simulation in section 4.1

``` r
library(TLOHO)
library(fields)
data = generate_simdata(n=100, rhox = 0, SNR = 4, option = "twoclusters")

fields::image.plot(matrix(data$beta.true, 30, 30)) # true beta image

graph0 = igraph::make_lattice(c(30,30))# construct lattice graph

fit <- tloho_lm(data$Y, data$X, graph0, Dahl = T) # fit

image.plot(matrix(fit$median_beta_est, 30, 30))
image.plot(matrix(fit$cluster_est_Dahl, 30, 30))
```

``` r
# normal means model when X = I
Y = data$beta.true + rnorm(900, sd = 0.5)

image.plot(matrix(Y, 30, 30)) # noisy data
fit2 <- tloho_lm_normalmeans(Y, graph0, Dahl = T)

image.plot(matrix(fit2$median_beta_est, 30, 30))
image.plot(matrix(fit2$cluster_est_Dahl, 30, 30))
```

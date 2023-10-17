
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TLOHO

<img src="TLoHo_longvertical_text.png" width="200">
<!-- badges: start --> <!-- badges: end -->

TLOHO is an accompanying R package of the paper:

Lee, C. J., Luo, Z. T., & Sang, H. (2021). T-LoHo: A Bayesian
Regularization Model for Structured Sparsity and Smoothness on Graphs.
*Advances in Neural Information Processing Systems 34 (NeurIPS 2021)*
[\[paper
link\]](https://proceedings.neurips.cc/paper/2021/hash/05a70454516ecd9194c293b0e415777f-Abstract.html)

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("changwoo-lee/TLOHO")
```

## Example

##### Running simulation in Section 4.1

``` r
library(TLOHO)
library(fields)
data = generate_simdata(n=100, rhox = 0, SNR = 4, option = "twoclusters", seed = 1) # change seed if needed

str(data$Y) # response is n = 100 vector
str(data$X) # predictor is n = 100 by p = 900 matrix
str(data$beta.true) # p=900 dimensional true coefficient (unknown)

fields::image.plot(matrix(data$beta.true, 30, 30)) # true beta image, assumed to be structured with many zeros

graph0 = igraph::make_lattice(c(30,30))# construct lattice graph

# fit the model
fit <- tloho_lm(data$Y, data$X, intercept = F, scale = T, graph0 = graph0, nsave = 1000, nburn = 40000, nthin = 10) 

image.plot(matrix(fit$median_beta_est, 30, 30))
image.plot(matrix(fit$cluster_est, 30, 30))
```

##### Running a model with an intercept term

``` r
library(TLOHO)
library(fields)

packageVersion("TLOHO") # make sure you are using version >1.2.0 
# generate non-scaled covariates
# here I set n=300, and p=900
data = generate_simdata(n=300, rhox = 0, SNR = 4, option = "twoclusters", seed = 1, scale = F) # change seed if needed
Y = data$Y + 2# shift response to have intercept 2
X = data$X
# not scaled
summary(colSums(X^2))

str(data$beta.true) # p=900 dimensional true coefficient (unknown)  WITHOUT intercept
fields::image.plot(matrix(data$beta.true, 30, 30)) # true beta image, assumed to be structured with many zeros WITHOUT intercept
graph0 = igraph::make_lattice(c(30,30))# construct lattice graph, WITHOUT intercept

# fit the model
# see ?tloho_lm for details on the function arguments
fit <- tloho_lm(Y, X, intercept = T, scale = F, graph0 = graph0, nsave = 1000, nburn = 40000, nthin = 10) 

# intercept
plot(fit$beta_out[,1], type="l")
abline(h=2, col = 2)

# coefficients (without intercept)
image.plot(matrix(fit$median_beta_est[-1], 30, 30))
image.plot(matrix(fit$cluster_est[-1], 30, 30))
```

##### Normal means model example (Section 4.2), on the lattice graph

``` r
# normal means model when X = I
Y = data$beta.true + rnorm(900, sd = 0.5)

image.plot(matrix(Y, 30, 30)) # noisy data
X = diag(900)

fit2 <- tloho_lm(Y, X, graph0 = graph0)

image.plot(matrix(fit2$median_beta_est, 30, 30))
image.plot(matrix(fit2$cluster_est, 30, 30))
```

##### Horseshoe+ prior instead of horseshoe (experimental)

``` r
# using data from the first example (Simulation data)
fit_hsplus <- tloho_lm(data$Y, data$X, graph0 = graph0, hsplus =T, seed = 123) # added in ver 1.1.0 (experimental)


image.plot(matrix(fit_hsplus$median_beta_est, 30, 30))


node = 1 # node that true beta = 0
plot(density(fit_hsplus$beta_out[,node]), main = "posterior of beta, black:hsplus, red:hs ")
lines(density(fit$beta_out[,node]), col = "red")
# hsplus is often more highly peaked at 0 compared to hs, but not always..
```

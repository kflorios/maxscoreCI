# maxscoreCI

`maxscoreCI` is an R package for maximum score estimation and confidence interval construction in binary response models. The package implements point estimation via differential evolution optimization (`DEoptim`) combined with C++ accelerated evaluation routines (`Rcpp`), and constructs valid 95% confidence intervals using the reshaped bootstrap method for estimators with cube-root asymptotics.

---

## Theoretical Background

Standard bootstrap methods often fail for non-smooth, cube-root asymptotic estimators such as Manski's maximum score estimator. This package implements the **Reshaped Bootstrap** methodology developed in:

> **Cattaneo, M. D., Jansson, M., & Nagasawa, K. (2020).**  
> *Bootstrap‐based inference for cube root asymptotics.*  
> **Econometrica**, 88(5), 2203–2219.  
> DOI: [10.3982/ECTA17950](https://doi.org/10.3982/ECTA17950)

The reshaped bootstrap corrects the asymptotic distribution by introducing a quadratic centering term driven by a smoothed sample Hessian matrix estimated via numerical derivatives or plug-in kernel estimation. The main result is offered with the plug-in kernel estimation method.

---

## Features

- **C++ Acceleration:** High-performance evaluation of maximum score objective functions and smoothed sample Hessians (`Rcpp`).
- **Global Optimization:** Differential Evolution optimization (`DEoptim`) to navigate discontinuous, non-smooth objective functions.
- **Normalization Control:** Option to specify the sign (`1` or `-1`) of the first covariate's coefficient for scale identification.
- **Reshaped Bootstrap Inference:** Automated point estimation and 95% confidence interval construction.

---

## Installation

You can install the package directly from GitHub using `devtools` or `remotes`:

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install maxscoreCI from GitHub
devtools::install_github("kflorios/maxscoreCI")
```

### Windows Prerequisites
Building C++ code on Windows requires **Rtools** installed matching your R version:
- Download and install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

---

## Quick Start Example

This example demonstrates how to perform estimation using sample data bundled with the package (`X.txt` and `y.txt`).

```r
library(maxscoreCI)
library(matrixStats)

# 1. Load sample dataset bundled with the package
x_path <- system.file("examples", "X.txt", package = "maxscoreCI")
y_path <- system.file("examples", "y.txt", package = "maxscoreCI")

X_raw <- read.table(x_path, header = FALSE)
y_raw <- read.table(y_path, header = FALSE)

# 2. Preprocess covariates and outcome
# Drop ID column (first column) and extract matrix X
X <- as.matrix(X_raw[, -1])

# Standardize covariates (except intercept/last column)
means <- colMeans(X)
sds   <- colSds(X)
for (j in 1:(ncol(X) - 1)) {
  X[, j] <- (X[, j] - means[j]) / sds[j]
}

# Extract binary response vector y
y <- as.numeric(y_raw[, 2])

# 3. Estimate Maximum Score & Compute 95% Confidence Intervals
# Note: Increase 'B' (e.g., B = 500) for empirical applications
set.seed(42)
results <- maxscore_ci(
  y = y, 
  X = X, 
  B = 8,            # Number of bootstrap replications
  plugin_bw = 2.0,  # Plug-in bandwidth choice
  num_eps = 0.5,    # Numerical derivative step size
  n_starts = 3,     # Global restarts for DEoptim
  x1_coef = 1       # Parameter fixing sign of X[, 1]
)

# 4. View Results
cat("\n=== Point Estimation and 95% Reshaped Bootstrap CIs ===\n")
print(results$confidence_interval_95)
```

---

## Function Reference

### `maxscore_ci(y, X, B = 50, plugin_bw = 2.0, num_eps = 0.5, n_starts = 3, x1_coef = 1)`

Main user-facing function to estimate parameters and construct confidence intervals.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `y` | Vector | *Required* | Binary dependent variable ($y \in \{0, 1\}$). |
| `X` | Matrix | *Required* | Covariate design matrix. |
| `B` | Integer | `50` | Number of bootstrap replications. |
| `plugin_bw` | Numeric | `2.0` | Bandwidth for kernel smoothing in plug-in Hessian calculation. |
| `num_eps` | Numeric | `0.5` | Epsilon step size for numerical derivative Hessian calculation. |
| `n_starts` | Integer | `3` | Number of random restarts for `DEoptim` search. |
| `x1_coef` | Numeric | `1` | Normalization parameter for the first covariate (`1` or `-1`). |

---



## License

This project is licensed under the MIT License.

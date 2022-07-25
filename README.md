# Layered_PCA
R implementation of Layered Principal Component Analysis (LPCA)

## Usage
```r
source("LPCA_majo.R")
X <- matrix(rnorm(1000), 100, 10)
res <- MULTIPLE_STARTS("LPCA_majo(X,2,3)",100)
```

## Reference
Yamashita, N. Principal component analysis constrained by layered simple structures. Adv Data Anal Classif (2022). https://doi.org/10.1007/s11634-022-00503-9

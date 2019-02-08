# glmnetPlus
glmnet package with user-provided initialization

Usage:

```r
library(glmnetPlus)

X <- matrix(rnorm(500*100), 500, 100)
y <- rnorm(500)

fit0 <- glmnetPlus::glmnet(X, y, thresh = 1e-20)
# initialize with the solution at lambda 50, and find solutions for lambda 51 to lambda 60
fitPlus <- glmnetPlus::glmnet(X, y, lambda = fit$lambda[51:60], beta0 = fit$beta[, 50], type.gaussian = "naive", thresh = 1e-20)
max(abs(fitPlus$beta - fit0$beta[, 51:60]))
```

best.lambda <- function(X1Y, ny, alpha) {
  lam.best <- vector()
  for (i in 1:1000) {
    train <- sample(1:nrow(X1Y), 30)
    lasso.md = glmnet(data.matrix(X1Y[-train,10:18]), as.numeric(X1Y[-train,ny]), alpha = alpha)
    pred = predict(lasso.md, data.matrix(X1Y[train,10:18]))
    rmse= sqrt(apply((X1Y[train,ny]-pred)^2,2,mean))
    lam.bes=lasso.md$lambda[order(rmse)[1]]
    lam.best[i] <- drop(lam.bes)
    lam.best[which(is.nan(lam.best))] <- 0
    }
  best <- mean(lam.best)
  if (best <= 0) {best = 0}
  return(best)
}
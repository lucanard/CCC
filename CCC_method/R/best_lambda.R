best_lambda <- function(X1Y, ny, alpha) {
  lam.best <- vector()
  for (i in 1:1000) {
    train <- sample(1:nrow(X1Y), 30)
    lasso.md = glmnet(X1Y[-train,10:18],X1Y[-train,ny], alpha = alpha)
    pred = predict(lasso.md,X1Y[train,10:18])
    rmse= sqrt(apply((X1Y[train,ny]-pred)^2,2,mean))
    lam.bes=lasso.md$lambda[order(rmse)[1]]
    lam.best[i] <- drop(lam.bes)
    #coef(lasso.md,s=lam.best)
  }
  best <- mean(lam.best)
  return(best)
}
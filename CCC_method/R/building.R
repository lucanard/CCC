building <- function(X1Y, ny, model) {
  if (model == "logistic") {
    Xa <- as.data.frame(X1Y[,10:18])
    YY <- as.numeric(X1Y[,ny])
    #if (max(YY) != 1) {stop("response must be logistic")}
    model <- glm(YY ~ .,  data = cbind.data.frame(YY, Xa), family = "binomial")
  }
  if (model == "pls") {
    pena <- penalized.pls.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = seq(-1000, 1000, 1), k=10, scale = T)
    lambda <- pena$lambda.opt
    ncomp <- pena$ncomp.opt
    model <- penalized.pls.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = lambda, ncomp = ncomp, k=10)
  }
  
  if (model == "lasso") {
    lambda <- best.lambda(X1Y, ny=ny, alpha = 1)
    model = glmnet(data.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), alpha = 1, lambda = lambda)
  }
  if (model == "ridge") {
    lambda <- best.lambda(X1Y, ny=ny, alpha = 0)
    model = glmnet(data.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), alpha = 0, lambda = lambda)
  }
  if (model == "b_ppls") {
    pena <- ppls.splines.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = seq(-1000, 1000, 1), k=10, scale = T, reduce.knots= TRUE)
    lambda <- pena$lambda.opt
    ncomp <- pena$ncomp.opt
    dummy <- X2s(as.matrix(X1Y[,10.18]), reduce.knots = TRUE)
    P <- Penalty.matrix(m = ncol(dummy$Z))
    model <- penalized.pls.cv(as.matrix(dummy$Z), P = P, lambda = lambda, as.numeric(X1Y[,ny]), k=10)
  }
  return(model)
}
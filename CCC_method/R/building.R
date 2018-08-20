building <- function(X1Y, ny, model) {
  if (model == "logistic") {
    Xa <- as.data.frame(X1Y[,10:18])
    YY <- as.numeric(X1Y[,ny])
    #if (max(YY) != 1) {stop("response must be logistic")}
    model <- glm(YY ~ .,  data = cbind.data.frame(YY, Xa), family = "binomial")
  } else if (model == "pls") {
    pena <- ppls::penalized.pls.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = seq(-1000, 1000, 1), k=10, scale = T)
    lambda <- pena$lambda.opt
    ncomp <- pena$ncomp.opt
    model <- ppls::penalized.pls.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = lambda, ncomp = ncomp, k=10)
  } else if (model == "lasso") {
    lambda <- best.lambda(X1Y, ny=ny, alpha = 1)
    model = glmnet::glmnet(data.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), alpha = 1, lambda = lambda)
  } else if (model == "ridge") {
    lambda <- best.lambda(X1Y, ny=ny, alpha = 0)
    model = glmnet::glmnet(data.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), alpha = 0, lambda = lambda)
  }else  if (model == "b_ppls") {
    pena <- ppls::ppls.splines.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = seq(-1000, 1000, 10), k=10, scale = T, reduce.knots= TRUE)
    lambda <- pena$lambda.opt
    ncomp <- pena$ncomp.opt
    dummy <- ppls::X2s(as.matrix(X1Y[,10:18]), reduce.knots = TRUE)
    P <- ppls::Penalty.matrix(m = ncol(dummy$Z))
    model <- ppls::penalized.pls.cv(as.matrix(dummy$Z), P = P, ncomp = ncomp, lambda = lambda, as.numeric(X1Y[,ny]), k=10)
  }
  return(model)
}
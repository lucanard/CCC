model.testing <- function(X1Y, ntest, ny, model, errori, per.test = FALSE) {
  train <- lapply(seq(1,ntest), function(x) sample(1:nrow(X1Y), 30))
  X <- lapply(train, function(x) X1Y[-x,10:18])
  X1 <- lapply(train, function(x) X1Y[x,10:18])
  Y <- lapply(train, function(x) X1Y[-x,ny])
  Ytest <- lapply(train, function(x) X1Y[x,ny])
  lambda = seq(-1000, 1000, 10)
  XX <- do.call(rbind, X1)
  ppms <- c(-50, -30, -10, -5, -3, 0, 3, 5, 10, 30, 50)
  Cees <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5)
  XXX <- list(length = errori)
  Ytest <- permutize(Ytest, per.test)
  X <- lapply(X, function(x) as.matrix(x))
  Y <- lapply(Y, function(x) as.numeric(x))
  Ytest <- lapply(Ytest, function(x) as.numeric(x))
  if (model == "lasso") {
    lambda <- best.lam(X1Y, ny=i, alpha = 1)}
  if (model == "ridge") {
    lambda <- best.lam(X1Y, ny=i, alpha = 0)}
  if (model == "pls") {
    pena <- penalized.pls.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = lambda, k=10, scale = T)
    lambda <- pena$lambda.opt
    ncomp <- pena$ncomp.opt
  }
  if (model == "b_ppls") {
    pena <- ppls.splines.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = lambda, k=10, scale = T, reduce.knots= TRUE)
    lambda <- pena$lambda.opt
    ncomp <- pena$ncomp.opt
  }
  rms <- vector(mode = "numeric", length = length(errori))
  means_allA <- vector(mode = "numeric", length = length(errori))
  sds_allA <- vector(mode = "numeric", length = length(errori))
  for (i2 in 1:length(ppms)) {
    XX1 <- errorize(XX, errori, i2)
    XX11 <- split(as.data.frame(XX1), rep(1:ntest, each=30))
    XXX <- lapply(XX11, function(x) as.matrix(x))
    rmses <- vector(mode = "numeric", length = ntest)
    resA <- vector(mode = "numeric", length = ntest)
    for (i in 1:ntest) {
      if (model == "plsr") {
        C <- data.frame(h = I(as.matrix(Y[[i]])), c = I(as.matrix(X[[i]])))
        PLSaaa <- plsr(method = c("oscorespls"), h ~ c, data = C, center = T, scale = T,
                       stdize = F, validation = "CV")
        PLStest1 <- data.frame(h = I(as.matrix(Ytest[[i]])), c = I(as.matrix(XXX[[i]])))
        test <- predict(PLSaaa, PLStest1, na.action = na.exclude)
        test <- round(as.data.frame(test))
      }
      if (model == "pls") {
        model.obje <- penalized.pls.cv(as.matrix(X[[i]]), as.numeric(Y[[i]]), lambda = lambda, ncomp = ncomp, k=10)
        testi <- new.penalized.pls(model.obje, as.matrix(XXX[[i]]))
        test <- testi$ypred}
      if (model == "b_ppls") {
        dummy <- X2s(as.matrix(X[[i]]), as.matrix(XXX[[i]]), reduce.knots = TRUE)
        P <- Penalty.matrix(m = ncol(dummy$Z))
        model.obje <- penalized.pls.cv(as.matrix(dummy$Z), P = P, lambda = lambda, as.numeric(Y[[i]]), k=10)
        testi <- new.penalized.pls(model.obje, as.matrix(dummy$Ztest))
        test <- testi$ypred}
      if (model == "lasso") {
        lasso.md = glmnet(as.matrix(X[[i]]), as.numeric(Y[[i]]), alpha = alpha, lambda = lambda)
        pred = predict(lasso.md,as.matrix(XXX[[i]]))
        test <- round(pred)
      }
      if (model == "ridge") {
        lasso.md = glmnet(as.matrix(X[[i]]), as.numeric(Y[[i]]), alpha = alpha, lambda = lambda)
        pred = predict(lasso.md,as.matrix(XXX[[i]]))
        test <- round(pred)
      }
      if (model == "logistic"){
        Xa <- as.data.frame(X[[i]])
        YY <- Y[[i]]
        if (max(as.numeric(Y[[i]])) != 1)
          stop("response must be logistic")
        bin.model <- glm(YY ~ .,  data = cbind.data.frame(YY, Xa), family = "binomial")
        Xa <- as.data.frame(XXX[[i]])
        pred = predict(bin.model, Xa, type="response")
        test <- round(pred)
      }
      test <- round(test)
      test[test < 0] <- 0
      test <- as.matrix(test)
      rmse= sqrt(apply((Ytest[[i]]-test)^2,2,mean))
      rmses[i] <- rmse
      result <- abs(test-Ytest[[i]])
      finres <- result == 0
      finres <- colSums(finres, "TRUE")
      resA[i] <- finres
    }
    rms[i2] <- mean(rmses)
    pper.resmeanA <- (mean(resA)/30)*100
    pper.ressdA <- (sd(resA)/30)*100
    means_allA[i2] <- pper.resmeanA
    sds_allA[i2] <- pper.ressdA
  }
  result <- cbind(means_allA, sds_allA, rms)
  return(result)
}
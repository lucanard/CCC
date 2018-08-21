model.test <- function(X1Y, ntest = ntest, model, errors, nyc) {
  nyc <- c(1,2,4,5,6,8,9)
  roc.perf <- list()
  roc.auc <- list()
  for (i3 in nyc){
    print(paste("beginning prediction of", colnames(X1Y)[i3]))
    ny = i3
    train <- lapply(seq(1,ntest), function(x) sample(1:nrow(X1Y), 30))
    X <- lapply(train, function(x) X1Y[-x,10:18])
    X1 <- lapply(train, function(x) X1Y[x,10:18])
    Y <- lapply(train, function(x) X1Y[-x,ny])
    Ytest <- lapply(train, function(x) X1Y[x,ny])
    lambda = seq(0, 1000, 10)
    XX <- do.call(rbind, X1)
    ppms <- c(-50, -30, -10, -5, -3, 0, 3, 5, 10, 30, 50)
    Cees <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5)
    if (errors == "ppms") {error <- ppms}
    if (errors == "Cees") {error <- Cees}
    if (errors == "nerr") {error <- vector(length=length(ppms))}
    if (errors == "none") {error <- vector(length = 1)}
    XXX <- list(length = error)
    X <- lapply(X, function(x) as.matrix(x))
    Y <- lapply(Y, function(x) as.numeric(x))
    Ytest <- lapply(Ytest, function(x) as.numeric(x))
    if (model == "lasso") {
      lambda <- best.lambda(X1Y, ny=ny, alpha = 1)
    }
    if (model == "ridge") {
      lambda <- best.lambda(X1Y, ny=ny, alpha = 0)
    }
    if (model == "pls") {
      pena <- ppls::penalized.pls.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = lambda, k=10, scale = T)
      lambda <- pena$lambda.opt
      ncomp <- pena$ncomp.opt
    }
    if (model == "b_ppls") {
      pena <- ppls::ppls.splines.cv(as.matrix(X1Y[,10:18]), as.numeric(X1Y[,ny]), lambda = lambda, k=10, scale = T, reduce.knots= TRUE)
      lambda <- pena$lambda.opt
      ncomp <- pena$ncomp.opt
    }
    resAA <- list()
    testY <- list()
    testY1 <- list()
    for (i2 in 1:length(error)) {
      if (errors != "none") {XX1 <- errorize(XX, errors, i2)} else {XX1 <- XX}
      XX11 <- split(as.data.frame(XX1), rep(1:ntest, each=30))
      XXX <- lapply(XX11, function(x) as.matrix(x))
      resA <- list(length = ntest)
      for (i in 1:ntest) {
        if (model == "pls") {
          model.obje <- ppls::penalized.pls.cv(as.matrix(X[[i]]), as.numeric(Y[[i]]), lambda = lambda, ncomp = ncomp, k=10)
          testi <- ppls::new.penalized.pls(model.obje, as.matrix(XXX[[i]]))
          test <- testi$ypred
        }
        if (model == "b_ppls") {
          dummy <- ppls::X2s(as.matrix(X[[i]]), as.matrix(XXX[[i]]), reduce.knots = TRUE)
          P <- ppls::Penalty.matrix(m = ncol(dummy$Z))
          model.obje <- ppls::penalized.pls.cv(as.matrix(dummy$Z), P = P, lambda = lambda, as.numeric(Y[[i]]), k=10)
          testi <- ppls::new.penalized.pls(model.obje, as.matrix(dummy$Ztest))
          test <- testi$ypred
        }
        if (model == "lasso") {
          lasso.md = glmnet::glmnet(data.matrix(X[[i]]), as.numeric(Y[[i]]), alpha = 1, lambda = lambda)
          pred = predict(lasso.md, data.matrix(XXX[[i]]))
          test <- pred
        }
        if (model == "ridge") {
          lasso.md = glmnet::glmnet(data.matrix(X[[i]]), as.numeric(Y[[i]]), alpha = 0, lambda = lambda)
          pred = predict(lasso.md, data.matrix(XXX[[i]]))
          test <- pred
        }
        if (model == "logistic"){
          Xa <- as.data.frame(X[[i]])
          YY <- Y[[i]]
          if (max(as.numeric(Y[[i]])) != 1)
            stop("response must be logistic")
          bin.model <- glm(YY ~ .,  data = cbind.data.frame(YY, Xa), family = "binomial")
          Xa <- as.data.frame(XXX[[i]])
          pred = predict(bin.model, Xa, type="response")
          test <- pred
        }
        #test <- round(test)
        test[test < 0] <- 0
        test <- as.matrix(test)
        resA[[i]] <- test
      }
      resAA[[i2]] <- t(do.call(cbind, resA))
      testY[[i2]] <- do.call(rbind, Ytest)
      Ytest1 <- permutize(Ytest)
      testY1[[i2]] <- do.call(rbind, Ytest1)
    }
    resAB <- do.call(rbind, resAA)
    testYB <- do.call(rbind, testY)
    manypred = ROCR::prediction(as.numeric(resAB), as.numeric(testYB))
    roc.perf[[i3]] = ROCR::performance(manypred, measure = "tpr", x.measure = "fpr")
    roc.auc[[i3]] = ROCR::performance(manypred, measure = "auc")
  }
  roc <- list(roc.perf, roc.auc) 
  names(roc) <- c("roc.perf", "roc.auc")
  roc1 <- roc
  
Sul_roc <- roc1$roc.perf[[1]]
pho_roc <- roc1$roc.perf[[2]]
acidic_roc <- roc1$roc.perf[[4]]
NN_roc <- roc1$roc.perf[[5]]
het_roc <- roc1$roc.perf[[6]]
CO_roc <- roc1$roc.perf[[8]]
bs_roc <- roc1$roc.perf[[9]]

sul_auc <- roc1$roc.auc[[1]]
pho_auc <- roc1$roc.auc[[2]]
acidic_auc <- roc1$roc.auc[[4]]
NN_auc <- roc1$roc.auc[[5]]
het_auc <- roc1$roc.auc[[6]]
CO_auc <- roc1$roc.auc[[8]]
bs_auc <- roc1$roc.auc[[9]]
auc <- c(sul_auc@y.values, pho_auc@y.values, 
              acidic_auc@y.values, NN_auc@y.values, 
              het_auc@y.values, CO_auc@y.values,
              bs_auc@y.values)
names(auc) <- c("sul", "pho", "acidic", "NN", "het", 
                     "CO", "bs")
results <- list(roc1, auc)

if (errors == "none") {sub_title = "with no measurement error"}
if (errors == "nerr") {sub_title = "with normal distributed error"}

color <- c("blue", "green", "orange", "dark red", "gray", "violet", "black", "red")
parameters <- c("sul", "pho", "acidic", "NN", "het", "CO", "bs", "random_pred")
plot(Sul_roc, col = "blue", main = paste("ROC curves", model, "model"), cex = 2, lwd = 2, lty = 1, sub = sub_title)
plot(pho_roc, add = T, col = "green", lwd = 2, lty = 2) 
plot(acidic_roc, add = T, col = "orange", lwd = 2, lty = 3)
plot(NN_roc, add = T, col = "dark red", lwd = 2, lty = 4)
plot(het_roc, add = T, col = "gray", lwd = 2, lty = 5)
plot(CO_roc, add = T, col = "violet", lwd = 2, lty = 6)
plot(bs_roc, add = T, col = "black", lwd = 2, lty = 3)
abline(a=0, b= 1, col = "red", lwd = 1)
legend("bottomright",legend=parameters, col=color, lwd=3, lty=c(1, 2, 3, 4, 5, 6, 3), cex=0.8, inset=c(0,0))
return (results)
}

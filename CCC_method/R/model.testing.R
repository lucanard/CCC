model.testing <- function(X1Y, ntest, nyc, model, errors) {
  options(warn=-1)
  results_phenolics <- list()
  results_aliph <- list()
  results <- list()
  for (i4 in nyc) {
    print(paste("beginning prediction of", colnames(X1Y)[i4]))
    ny <- i4  
    train <- lapply(seq(1,ntest), function(x) sample(1:nrow(X1Y), 30))
    X <- lapply(train, function(x) X1Y[-x,10:18])
    X1 <- lapply(train, function(x) X1Y[x,10:18])
    Y <- lapply(train, function(x) X1Y[-x,ny])
    Ytest <- lapply(train, function(x) X1Y[x,ny])
    lambda = seq(-1000, 1000, 10)
    XX <- do.call(rbind, X1)
    ppm <- c(-50, -30, -10, -5, -3, 0, 3, 5, 10, 30, 50)
    Cees <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5)
    if (errors != "ppm" & errors != "Cees" & errors != "nerr" & errors != "none") {errors = "none"}
    if (errors == "ppm") {error <- ppm}
    if (errors == "Cees") {error <- Cees}
    if (errors == "nerr") {error <- vector(length=length(ppm))}
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
    nami <- vector(length = length(errors))
    MSE_permuted <- vector(mode = "numeric", length = length(errors))
    TPR_Permuted <- vector(mode = "numeric", length = length(errors))
    TPR_sd_permuted <- vector(mode = "numeric", length = length(errors))
    MSE <- vector(mode = "numeric", length = length(errors))
    TPR <- vector(mode = "numeric", length = length(errors))
    TPR_sd <- vector(mode = "numeric", length = length(errors))
    for (i2 in 1:length(error)) {
      if (errors != "none") {XX1 <- errorize(XX, errori = errors, i2)} else {XX1 <- XX}
      XX11 <- split(as.data.frame(XX1), rep(1:ntest, each=30))
      XXX <- lapply(XX11, function(x) as.matrix(x))
      MSEs <- vector(mode = "numeric", length = ntest)
      resA <- vector(mode = "numeric", length = ntest)
      MSEs1 <- vector(mode = "numeric", length = ntest)
      resA1 <- vector(mode = "numeric", length = ntest)
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
          test <- round(pred)
        }
        if (model == "ridge") {
          lasso.md = glmnet::glmnet(data.matrix(X[[i]]), as.numeric(Y[[i]]), alpha = 0, lambda = lambda)
          pred = predict(lasso.md, data.matrix(XXX[[i]]))
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
        if (model == "logistic") {
          MSEi = mean((abs(Ytest[[i]]-test))^2)
         } else {
        MSEi = apply((abs(Ytest[[i]]-test))^2,2,mean)
         }
        MSEs[i] <- MSEi
        Ytest1 <- permutize(Ytest)
        if (model == "logistic") {
          MSE1 = mean((abs(Ytest1[[i]]-test)^2))
        } else {
        MSE1 = apply((abs(Ytest1[[i]]-test)^2),2,mean)
        }
        MSEs1[i] <- MSE1
        test <- round(test)
        test[test < 0] <- 0
        test <- as.matrix(test)
        result <- abs(test-Ytest[[i]])
        finres <- result == 0
        finres <- colSums(finres, "TRUE")
        resA[i] <- finres
        result1 <- abs(test-Ytest1[[i]])
        finres1 <- result1 == 0
        finres1 <- colSums(finres1, "TRUE")
        resA1[i] <- finres1
      }
      MSE[i2] <- mean(MSEs)
      pper.resmeanA <- (mean(resA)/30)*100
      pper.ressdA <- (sd(resA)/30)*100
      TPR[i2] <- pper.resmeanA
      TPR_sd[i2] <- pper.ressdA
      MSE_permuted[i2] <- mean(MSEs1)
      pper.resmeanA <- (mean(resA1)/30)*100
      pper.ressdA <- (sd(resA1)/30)*100
      TPR_Permuted[i2] <- pper.resmeanA
      TPR_sd_permuted[i2] <- pper.ressdA
      nami[i2] <- as.character(paste0(errors, "_", error[i2]))
    }
      tutti <- cbind(TPR, TPR_sd, MSE, TPR_Permuted, TPR_sd_permuted, MSE_permuted)
      row.names(tutti) <- nami
      results[[i4]] <- tutti
    if (errors == "ppm" | errors == "Cees") {
   if (errors == "ppm") {sub_title = "with mass measurement error"}
    if (errors == "Cees") {sub_title = "with Carbon estimation errors"}
      err = error
    true_positive_rate <- TPR
      perm.test <- TPR_Permuted
      plot(err, true_positive_rate, type="l", col = "black", ylim=c(min((sapply(perm.test, function(x) min(x)))),100), main = paste("True positive rate curve of", names(X1Y[1,][ny])), sub = sub_title, xlab = errors, lwd = 2)
      lines(err, true_positive_rate, lty = 3, col = "black")
      lines(err, true_positive_rate + TPR_sd,lty=2, col = "grey")
      lines(err, true_positive_rate - TPR_sd,lty=2, col = "grey")
      lines(err, perm.test, type = "l", col = "red", lwd = 2)
      abline(v=0, h=NULL, col="blue", lwd = 2)
      abline(v=-(max(err)*(1/3)), h=NULL, col="green")
      abline(v=(max(err)*(1/3)), h=NULL, col="green")
      legend("right",legend=c("TPR", "permuted TPR", "TPR_SD"), col= c("black", "red", "grey"), lwd=2, cex=0.6)
  }
  }
  result = results[-which(sapply(results, is.null))]
  return(result)
}
  

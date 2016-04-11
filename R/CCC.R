CCC_code <- function(x) {
  mat <- function(x){
    phe <- character(length=length(x))
    for (i in 1:length(x))
      if (x[[i]][1][[1]] == TRUE)
      {phe[i] <- length(x[[i]][2][[1]])}
    else
    {phe[i] <- 0 }
    phe}
  sulfur <- str_detect(x, pattern="[Ss]") + 0
  Carbon <- str_count(x, pattern="[Cc]") + 0
  Oxygen <- str_count(x, pattern="[Oo]") + 0
  CO <- as.numeric((Carbon)/(Oxygen))
  CO[CO < 7] <- 0
  CO[CO >= 7] <- 1
  NN <- str_detect(x, pattern="[Nn]")
  phosphate <- as.numeric(str_detect(x, pattern = "[Pp]"))
  CbC <- str_count(x, fixed("CC=CC"))
  CaC <- str_count(x, fixed("C\\C=C/C"))
  CC <- str_count(x, fixed("CC(C)C"))
  CCC <- str_count(x, fixed("CCCC"))
  aliph <- as.numeric(CbC + CaC + CC + CCC)
  A <- smiles2sdf(x)
  Pyri = "c1ccncc1"
  Pyrr = "c1cc[nH]c1"
  fur = "c1ccoc1"
  ncyc = "c1cncn1"
  query <- 'c1ccccc1'
  query1 <- "OCC(O)CO"
  query2 <- "O=CO"
  sdfset <- smiles2sdf(c(Pyri, Pyrr, fur, ncyc))
  sdfset@ID <- c("Pyri", "Pyrr", "fur", "ncyc")
  mcs <- fmcsBatch(sdfset, A, al = 0, au = 0, bl= 0, bu = 0, matching.mode = "aromatic", numParallel = 2)
  mcs[mcs[,5] < 1] <- 0
  het <- drop(mcs[,5])
  mols <- tryCatch({sapply(x, parse.smiles)}, error = function (e) {sapply(x, function(x) parse.smiles(x, kekulise = FALSE))})
  does2 <- matches(query2, mols, return.matches = FALSE)
  does <- matches(query, mols, return.matches = TRUE)
  does1 <- matches(query1, mols, return.matches = FALSE)
  acidic <- as.numeric(does2 + 0)
  bs <- as.numeric(does1 + 0)
  phenolics <- as.numeric(mat(does))
  Y <- cbind(sulfur, phenolics, acidic, NN, aliph, CO, bs)
  return(Y)
}
reorder.csv <- function(directory, compi, selection) {
    setwd(directory)
  file <- list.files(pattern = "csv", full.names = FALSE)
  candidate <- list()
  for (i in 1:length(file)) {
    print(file[i])
    files <- read.csv(file[i], stringsAsFactors=FALSE)
    name <- files[4,]
    files <- files[5:nrow(files),]
    names(files) <- name
    es <- files$Smiles
    com <- CCC_code(es)
    predicted <- compi[selection, 10:16]
    ebb <- apply(com, 1, '-', as.numeric(predicted))
    eb <- colSums(abs(ebb))
    eb <- eb + 1
    reordi <- order(eb)
    correct <- files[reordi,]
    candidate[[i]] <- correct
    names(candidate[[i]]) <- names(files)
  }
  return(candidate)
}
reorder.sdf <- function(directory, compi, selection) {
  setwd(directory)
  file <- list.files(pattern = "sdf", full.names = FALSE)
  fil <- str_split(file, pattern = "_")
  fi <- do.call(cbind, fil)
  fil <- fi[1,]
  candidate <- list()
  for (i in 1:length(file)) {
    print(fil[i])
    esp <- read.SDFstr(file[i])
    esp <- as(esp, "SDFset")
    espo <- sdf2smiles(esp)
    es <- toString(espo)
    es <- as.vector(unlist(strsplit(es, split=", ", perl = TRUE)))
    com <- CCC_code(es)
    predicted <- compi[selection, 10:16]
    ebb <- apply(com, 1, '-', as.numeric(predicted))
    eb <- colSums(abs(ebb))
    eb <- eb + 1
    reordi <- order(eb)
    correct <- esp[reordi,]
    candidate[[i]] <- correct
  }
  return(candidate)
}
best_lambda <- function(X1Y, ny, alpha) {
  require(glmnet, quietly =T)
  lam.best <- vector()
  for (i in 1:1000) {
    train <- sample(1:nrow(X1Y), 30)
    lasso.md = glmnet(X1Y[-train,10:18],X1Y[-train,ny], alpha = alpha)
    pred = predict(lasso.md,X1Y[train,10:18])
    rmse= sqrt(apply((X1Y[train,ny]-pred)^2,2,mean))
    lam.bes=lasso.md$lambda[order(rmse)[1]]
    lam.best[i] <- drop(lam.bes)
    rms <- rmse[order(rmse)[1]]
    rmses[i] <- drop(rms)
    #coef(lasso.md,s=lam.best)
  }
  best <- mean(lam.best)
  return(best)
}
dataset.building <- function(x) {
  require(stringr, quietly = T)
  require(ChemmineR, quietly = T)
  require(rcdk, quietly = T)
  require(fmcsR, quietly = T)
  is.odd <- function(x) x %% 2 != 0
  Pyri = "c1ccncc1"
  Pyrr = "c1cc[nH]c1"
  fur = "c1ccoc1"
  ncyc = "c1cncn1"
  query <- 'c1ccccc1'
  query1 <- "OCC(O)CO"
  query2 <- "O=CO"
  sdfset <- smiles2sdf(c(Pyri, Pyrr, fur, ncyc))
  sdfset@ID <- c("Pyri", "Pyrr", "fur", "ncyc")
  RT <- as.numeric(x[["RT"]])
  mass <- as.numeric(x[["MM"]])
  comp <- str_replace(x[["Formula"]], pattern ="C", replacement = "")
  com <- str_split(comp, pattern = "H", n = Inf)
  nC <- sapply(com, function(x) drop(as.numeric(x[1])))
  nC[is.na(nC)] <- 1
  nm <- round(mass)
  md <- (mass-nm)
  RMD <- (md/nm*1000000)
  mC <- (as.numeric(nC*12))
  pC <- (mC/(mass)*100)
  rmass <- mass-mC
  rnm <- round(rmass)
  rmd <- rmass-rnm
  rRMD <- rmd/rnm*1000000
  odd <- is.odd(nm)
  odd <- odd + 0
  Sulfur <- str_detect(x[["Formula"]], pattern = "S")+0
  mat <- function(x){
    phe <- character(length=length(x))
    for (i in 1:length(x))
      if (x[[i]][1][[1]] == TRUE)
      {phe[i] <- length(x[[i]][2][[1]])}
    else
    {phe[i] <- 0 }
    phe}
  phenolics <- vector(length=nrow(x),mode="numeric")
  het <- vector(length=nrow(x),mode="numeric")
  acidic <- vector(length=nrow(x),mode="numeric")
  bs <- vector(length=nrow(x),mode="numeric")
  CbC <- str_count(x$smiles, fixed("CC=CC"))
  CaC <- str_count(x$smiles, fixed("C\\C=C/C"))
  CC <- str_count(x$smiles, fixed("CC(C)C"))
  CCC <- str_count(x$smiles, fixed("CCCC"))
  aliph <- as.numeric(CbC + CaC + CC + CCC)
  A <- smiles2sdf(x$smiles)
  A@ID <- x$compound
  mcs <- fmcsBatch(sdfset, A, al = 0, au = 0, bl= 0, bu = 0, matching.mode = "aromatic", numParallel = 2)
  mcs[mcs[,5] < 1] <- 0
  het <- drop(mcs[,5])
  mols <- sapply(x$smiles, parse.smiles)
  does2 <- matches(query2, mols, return.matches = FALSE)
  acidic <- as.numeric(does2 + 0)
  does <- matches(query, mols, return.matches = TRUE)
  does1 <- matches(query1, mols, return.matches = FALSE)
  bs <- as.numeric(does1 + 0)
  phenolics <- as.numeric(mat(does))
  sulfur <- str_detect(x$smiles, pattern="[Ss]") + 0
  Carbon <- str_count(x$smiles, pattern="[Cc]") + 0
  Oxygen <- str_count(x$smiles, pattern="[Oo]") + 0
  CO <- as.numeric((Carbon)/(Oxygen))
  CO[CO < 7] <- 0
  CO[CO >= 7] <- 1
  NN <- str_detect(x$Formula, pattern = "N")
  phosphate <- as.numeric(str_detect(x$Formula, pattern = "P"))
  X1Y <- cbind(sulfur, phosphate, phenolics, acidic, NN, het, aliph, CO, bs, RT, mass, nC, md, RMD, pC, rRMD, odd, Sulfur)
  rownames(X1Y) <- x[["compound"]]
  return(X1Y)
}
errorize <- function (x, errori, i2) {
  X2 <- data.frame()
  Xm <- x
  if (errori == "ppms") {
    x[,2] <- (x[,2]+(x[,2]/1000000*ppms[i2]))
    x[,4] <- (x[,2] - round(x[,2]))
    x[,5] <- ((x[,4]/x[,2])*1000000)
    x[,7] <- ((x[,2]-(x[,3]*12) - round(x[,2]-(x[,3]*12)))/round(x[,2]-(x[,3]*12)))*1000000
  }
  if (errori == "Cees") {
    x[,3] <- (x[,3] + Cees[i2])
    maxC <- round(x[,2]/12)
    wC <- which(x[,3] - maxC >= 0)
    x[wC,3] <- maxC[wC] - 1
    x[which(x[,3] <= 0), 3] <- 1
    x[,6] <- (x[,3]*12/x[,2])*100
    x[,7] <- ((x[,2]-(x[,3]*12) - round(x[,2]-(x[,3]*12)))/round(x[,2]-(x[,3]*12)))*1000000
  }
  if (errori == "nerr") {
    x[,2] <- (x[,2]+rnorm(x[,2], mean=0, sd=0.001))
    x[,4] <- (x[,2] - round(x[,2]))
    x[,5] <- ((x[,4]/x[,2])*1000000)
    x[,3] <- round(x[,3]+rnorm(x[,3], mean=0, sd=1))
    maxC <- round(x[,2]/12)
    wC <- which(x[,3] - maxC >= 0)
    x[wC,3] <- maxC[wC] - 1
    x[which(x[,3] <= 0), 3] <- 1
    x[,6] <- ((x[,3]*12/x[,2])*100)
    x[,7] <- ((x[,2]-(x[,3]*12) - round(x[,2]-(x[,3]*12)))/(round(x[,2]-(x[,3]*12))))*1000000
  }
  if (errori == "none") {x <- x}
  print(sqrt(mean((x[,3]-Xm[,3])^2)))
  print(sqrt(mean((x[,2]-Xm[,2])^2)))
  print(sqrt(mean((x[,7]-Xm[,7])^2)))
  X2 <- x
  return(X2)
}
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
  if (model == "lasso") {require(glmnet, quietly = T)
    lambda <- best.lam(X1Y, ny=i, alpha = 1)}
  if (model == "ridge") {require(glmnet, quietly = T)
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
permutize <- function (x, per.test = FALSE){
  if (per.test != FALSE){
    x <- lapply(x, function(x) x[sample(length(x))])
  } else {x <- x}
  return(x)
}
apply.model <- function(xsaFA) {
  getIsoGroups=function(xsaFA){
    iso_groups=list()
    for (i in 1:length(xsaFA@isotopes)){
      temp=xsaFA@isotopes[[i]]$y
      if (!is.null(temp)){
        if (length(iso_groups)>=temp){
          iso_groups[[temp]]=c(iso_groups[[temp]],i)
        }else{
          iso_groups[[temp]]=i
        }
      }
    }
    return(iso_groups)
  }
  IsoGroup2ratios=function(xsaFA,isogroups,group_nr,peaklist){
    if (xsaFA@polarity=="positive"){mode=1}
    if (xsaFA@polarity=="negative"){mode=2}
    masses=xsaFA@groupInfo[isogroups[[group_nr]],'mz']
    camgroup <- vector()
    gro <- vector()
    data = peaklist[isogroups[[group_nr]],9:(ncol(peaklist)-3)]
    int = rowSums(peaklist[isogroups[[group_nr]], 9:(ncol(peaklist)-3)])
    gro = peaklist[isogroups[[group_nr]], ncol(peaklist)]
    camgroup = peaklist
    order=order(masses)
    masses=masses[order]+c(-1.007276,1.007276)[mode]
    data=data[order,]
    gro = gro[order]
    int = int[order]
    ratios=numeric()
    ratios[1]=1
    for (i in 2:length(masses)){
      to_rem = !(        ((data[1,]>0) & !is.na(data[1,]))              &        ((data[i,]>0)  & !is.na(data[i,]))       )
      m=data[1,!to_rem]
      m_iso=data[i,!to_rem]
      
      ratios[i] = sum(m_iso)/sum(m)
    }
    result=c()
    result$mz=masses
    result$ratios=ratios
    result$gro = gro
    result$int = int
    return(result)
  }
  rtsing <- function(xsaFA, intval = "into"){
    isogroups <- getIsoGroups(xsaFA)
    rts <- sapply(isogroups,function(x) mean(peaklist[x,"rt"]))
    return(rts)
  }
  rating <- function (xsaFA, intval = "into") {
    isogroups <- getIsoGroups(xsaFA)
    rts <- sapply(isogroups,function(x) mean(peaklist[x,"rt"]))
    ratios <- list()
    for(i in 1:length(isogroups)){
      ratios[[i]] <- IsoGroup2ratios(xsaFA,isogroups,group_nr = i, peaklist=peaklist)
    }
    return(ratios)
  }
  moverz <- function(ratios){
    mzs <- vector(length=length(ratios))
    for (i in 1:length(ratios)){
      mzs[[i]] <- round(tn$mz[[i]][1], digits=3)
    }
    return(mzs)
  }
  times <- function(ratios){
    RTS <- vector(length= length(ratios))
    for (i in 1:length(ratios)){
      RTS[[i]] <- drop(tn$V1[[i]][1])
    }
    return(RTS)
  }
  Car.est <- function(tn) {
    data(lin.numC)
    temp3 =vector()
    for (i in 1:nrow(tn)){
      temp3[[i]] <- as.vector(tn$ratios[[i]][2])/as.vector(tn$ratios[[i]][1])
    }
    C.est <- round(predict(lin.numC, newdata= data.frame(temp3),
                           interval = c("confidence"), level = 0.95))
    CC <- C.est[,1]
    return(CC)
  }
  querying <- function(tn, intval = "into"){
    structuring <- function(tn){
      RT <- as.numeric(tn$RTs)
      mass <- as.numeric(tn$mzs)
      nC <- as.numeric(tn$estC)
      nm <- round(mass)
      md <- (mass-nm)
      RMD <- (md/nm*1000000)
      mC <- (as.numeric(nC*12))
      pC <- (mC/(mass)*100)
      rmass <- mass-mC
      rnm <- round(rmass)
      rmd <- rmass-rnm
      rRMD <- rmd/rnm*1000000
      is.odd <- function(x) x %% 2 != 0
      odd <- is.odd(nm)
      odd <- odd + 0
      IMDP <- function(tn){
        IMDPs <- vector(length=nrow(tn))
        for (i in 1:nrow(tn)){if (length(tn$mz[[i]]) < 3) {IMDPs[[i]] <- 0}
          else {if (tn$mz[[i]][3]-(tn$mz[[i]][2]+1) < 0 &
                    tn$ratios[[i]][3]/tn$ratios[[i]][2] >= 0.26
                    & tn$ratios[[i]][3]/tn$ratios[[i]][2]<=0.75) {IMDPs[[i]] <- 1}
            else {IMDPs[[i]] <- 0}}}
        return(IMDPs)
      }
      IMD <- IMDP(tn)
      Sulfur <- IMD
      pX <- data.frame(RT, mass, nC, md, RMD, pC, rRMD, odd, Sulfur)
      data(bin.model.SS)
      data(pls.md.SS.Rdata)
      data(pls.md.phenolics)
      data(bin.model.acid.Rdata)
      data(bin.model.NN.Rdata)
      data(lasso.md.aliph.Rdata)
      data(bin.model.CO.Rdata)
      data(lasso.md.CO.Rdata)
      data(bin.model.bs.Rdata)
      SS <- predict(bin.model.SS, pX, type="response")
      acid <- predict(bin.model.acid, pX, type="response")
      NN <- predict(bin.model.NN, pX, type="response")
      CO <- predict(bin.model.CO, pX, type="response")
      bs <- predict(bin.model.bs, pX, type="response")
      phenolic <- new.penalized.pls(pls.md.phenolics, as.matrix(pX))
      phenolics <- phenolic$ypred
      names(phenolics) <- phenolics
      aliph <- predict(lasso.md.aliph, as.matrix(pX))
      strtn <- cbind(SS, phenolics, acid, NN, aliph, CO, bs)
      names(strtn) <- c("SS", "phenolics", "acid", "NN", "aliph", "CO", "bs")
      strtn <- as.data.frame(strtn)
      strtn <- round(strtn)
      strtn[strtn<0] <- 0
      names(strtn[2]) <- "phenolics"
      names(strtn[5]) <- "aliph"
      return(strtn)
    }
    grouping <- function (tn){
      gro <- lapply(tn$gro, function(x) drop(as.numeric(x[[1]])))
      gro <- as.numeric(gro)
      tn$gro <- gro
      return(tn)
    }
    ps_spec <- function(tn, peaklist){
      psspectra <- list()
      for (i in 1:length(tn$gro)) {
        mzo = peaklist$mz[which(peaklist$pcgroup == tn$gro [i])]
        int = rowSums(peaklist[which(peaklist$pcgroup == tn$gro [i]),9:(ncol(peaklist)-3)])
        temp1 = data.frame(mzo, int)
        psspectra[[i]] = drop(temp1)
      }
      names(psspectra) <- tn$gro
      order <- order(names(psspectra))
      psspectra <- unique(psspectra[order])
      return(psspectra)
    }
    
    strtn <- structuring(tn)
    tn <- grouping(tn)
    ps_spectra <- ps_spec(tn, peaklist)
    str <- do.call(list, strtn)
    tni <- cbind(tn, strtn)
    tni <- split(tni, f = tni$gro)
    return (tni)
  }
  if (str_detect(xsaFA@xcmsSet@filepaths[[1]], coll("neg")) == TRUE) {
    xsaFA@polarity <- "negative"} else {xsaFA@polarity <- "positive"}
  peaklist <- getPeaklist(xsaFA, intval='into')
  isogroups <- getIsoGroups(xsaFA)
  rts <- rtsing(xsaFA, intval = "into")
  ratios <- rating(xsaFA, intval="into")
  tn <- mapply(c, (rts/60), as.matrix(ratios), SIMPLIFY = "FALSE")
  tn <- as.data.frame(t(tn))
  mzs <- moverz(ratios)
  RTS <- times(ratios)
  CC <- Car.est(tn)
  tn <- cbind(mzs, tn, RTS, CC)
  colnames(tn)[1] <- "mzs"
  colnames(tn)[7] <- "RTs"
  colnames(tn)[8] <- "estC"
  colnames(tn)[2] <- "allRT"
  tni <- querying(tn)
  compi <- do.call(rbind, tni)
  return(compi)
}
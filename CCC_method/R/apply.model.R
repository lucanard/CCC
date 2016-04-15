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
#' @title Application of the CCC method on the xsaFA from CAMERA
#' @description the function uses the xsaFA given by CAMERA to group the isotopic patterns and it uses this data to give an estimation of their amount of Carbons and predicts the CCC Y dependent variables on each of the isotopic pattern. It can also handle statistical models built on different standards dataset using the function model.building.
#' @param xsaFA the xsaFA object resulting from CAMERA analysis 
#' @param models the statistical models to use in the CCC_method. Default is NULL, so the original models are loaded and applied. Elsewhere, own models can be applied , building them using the function "model.building".   
#' @return a list of grouped features with their isotopic pattern, estimated amount of Carbon and the CCC method predictions.
#' @usage CCC.xsaFA(xsaFA, models = NULL)
#' @export "CCC.xsaFA"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples 
#' #do not run: try only with your own xsaFA object, it requires xcms and CAMERA
#' tni <- CCC.xsaFA(xsaFA)
CCC.xsaFA <- function(xsaFA, models = NULL) {
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
  rowname <- row.names(peaklist[isogroups[[group_nr]], ])
  rts = mean(peakTable[isogroups[[group_nr]], "rt"])
  order=order(masses)
  masses=masses[order]+c(-1.007276,1.007276)[mode]
  data=data[order,]
  gro = gro[order]
  int = int[order]
  rowname = rowname[order]
  rts <- rts[order]
    ratios=numeric()
  ratios[1]=1
     for (i in 2:length(masses)){
    to_rem = !(        ((data[1,]>0) & !is.na(data[1,]))              &        ((data[i,]>0)  & !is.na(data[i,]))       )
    m=data[1,!to_rem]
    m_iso=data[i,!to_rem]

    ratios[i] = sum(m_iso)/sum(m)
  }
  result=c()
  result$rts=rts
  result$mz=masses
  result$ratios=ratios
  result$gro = gro
  result$int = int
  result$rowname <- rowname
  return(result)
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
Rt.est <- function(tn, peaklist) {
  if (str_detect(tail(colnames(peaklist), n = 1L), coll("RP")) & (max(peaklist[,"rt"]) <= 30)) {
    load(system.file("extdata", "rts.lm.rda", package = "CCC"))
    short <- as.numeric(sapply(tn$rts, "[[", 1))
    RT.est <- predict(rts.lm, newdata = data.frame(short), interval = c("confidence"), level = 0.95)
    RTS <- RT.est[,1]} else {RTS <- as.numeric(sapply(tn$rts, "[[", 1))}
  RTS[RTS>=60] <- 60
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
    if (is.null(models)) {
      load(system.file("extdata", "bin.model.acid.rda", package = "CCC"))
      load(system.file("extdata","bin.model.NN.rda", package ="CCC"))
      load(system.file("extdata","bin.model.SS.rda", package ="CCC"))
      load(system.file("extdata","bin.model.bs.rda", package ="CCC"))
      load(system.file("extdata","bin.model.CO.rda", package ="CCC"))
      load(system.file("extdata","lasso.md.CO.rda", package ="CCC"))
      load(system.file("extdata","lasso.md.aliph.rda", package ="CCC"))
      load(system.file("extdata","pls.md.SS.rda", package ="CCC"))
      load(system.file("extdata","pls.md.phenolics.rda", package ="CCC"))
    } else {
      bin.model.SS <- models[[1]]
      bin.model.acid <- models[[4]]
      bin.model.NN <- models[[5]]
      bin.model.CO <- models[[8]]
      bin.model.bs <- models[[9]]
      lasso.md.aliph <- models[[7]]
      pls.md.phenolics <- models[[3]]
    }
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
ratios <- rating(xsaFA, intval="into")
tn <- mapply(c, as.matrix(ratios), SIMPLIFY = "FALSE")
tn <- as.data.frame(t(tn))
tn$rowname <- lapply(tn$rowname, "[[", 1)
tn$rts <- sapply(tn$rts, "[[", 1)
mzs <- round(as.numeric(sapply(tn$mz, "[[", 1)), digits = 3)
RTS <- Rt.est(tn, peaklist)
CC <- Car.est(tn)
tn <- cbind(mzs, tn, RTS, CC)
colnames(tn)[1] <- "mzs"
colnames(tn)[8] <- "RTs"
colnames(tn)[9] <- "estC"
colnames(tn)[2] <- "allRT"
tni <- querying(tn)
return(tni)
}
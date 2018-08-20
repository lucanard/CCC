#' @title manual application of the CCC method on data
#' @description The function uses the users data to predict the CCC Y dependent variables
#' @param query a matrix data containing the RT, the mz of the main ions and their intensity
#' @param datapath the path to upload a csv files containing the RT, the mz of the main ions and their intensity 
#' @param polarity the polarity of the experiment
#' @param models the statisticals models to apply for the prediction. If NULL, old models from previous CCC version will be uploaded
#'
#' @return a data.frame containing the predicted value for your specific query 
#' @export "offline.CCC"
#' @author Luca Narduzzi "nardluca@gmail.com"
#'
#' @examples
#' load(system.file("extdata", "models.rda", package = "CCC"))
#' load(system.file("extdata", "example1.rda", package = "CCC"))
#' offline.CCC(query = example1, polarity = "negative", models = models)
offline.CCC <- function(query = NULL, datapath, polarity, models = NULL) {
  if (is.null(query)) {query <- read.csv(datapath, header = T, stringsAsFactors = F)} else {query = query}
  Car.est <- function(query) {
    load(system.file("extdata", "lin.numC.rda", package = "CCC"))
    temp3 <- query[2,3]/query[1,3]
    C.est <- round(predict(lin.numC, newdata= data.frame(temp3),
                           interval = c("confidence"), level = 0.95))
    CC <- C.est[,1]
    return(CC)
  }
  RT <- as.numeric(query[1,1])
  mass <- as.numeric(max(query$mz[which.max(query$int)]))
  if (polarity == "negative") {mass = mass + 1.0078}
  if (polarity == "positive") {mass = mass - 1.0078}
  if (is.numeric(polarity) == TRUE) {mass = mass - polarity}
  if (exists("polarity") == FALSE) {stop("please set experimental polarity")}
  nC <- as.numeric(Car.est(query))
  nm <- round(mass)
  maxC <- round(nm/12)
  wC <- which(nC - maxC >= 0)
  nC[wC] <- maxC[wC] - 1
  md <- (mass-nm)
  RMD <- (md/nm*1000000)
  mC <- (as.numeric(nC*12))
  pC <- (mC/(mass)*100)
  rmass <- mass-mC
  rnm <- round(rmass)
  rmd <- rmass-rnm
  rRMD <- rmd/rnm*1000000
  is.odd <- function(x) x %% 2 != 0
  odd <- is.odd(round(mass))
  odd <- odd + 0
  ratio <- query[3,3]/query[2,3]
  IMDPs <- function(tn){
    IMDPss <- vector(length=nrow(query))
    if (length(query$mz) < 3) {IMDPss <- 0}
      else {if (query$mz[3]-(query$mz[2]+1) < 0 &
                ratio >= 0.26
                & ratio <=0.75) {IMDPss <- 1}
        else {IMDPss <- 0}}
    return(IMDPss)
  }
  IMD <- IMDPs(tn)
  Sulfur <- IMD
  pX <- data.frame(RT, mass, nC, md, RMD, pC, rRMD, odd, Sulfur)
  if (is.null(models)) {
    load(system.file("extdata", "bin.model.acid.rda", package = "CCC"))
    load(system.file("extdata","b_pPLS.md.CO.rda", package ="CCC"))
    load(system.file("extdata","b_pPLS.md.aliph.rda", package ="CCC"))
    load(system.file("extdata","bin.model.NN.rda", package ="CCC"))
    load(system.file("extdata","bin.model.SS.rda", package ="CCC"))
    load(system.file("extdata","bin.model.bs.rda", package ="CCC"))
    load(system.file("extdata","bin.model.CO.rda", package ="CCC"))
    load(system.file("extdata","lasso.md.CO.rda", package ="CCC"))
    load(system.file("extdata","lasso.md.aliph.rda", package ="CCC"))
    load(system.file("extdata","pls.md.SS.rda", package ="CCC"))
    load(system.file("extdata","pls.md.phenolics.rda", package ="CCC"))
    models <- list(bin.model.SS, pls.md.phenolics, bin.model.NN, lasso.md.aliph, bin.model.CO, bin.model.bs)
  }else {
    bin.model.SS <- models[[1]]
    bin.model.pho <- models[[2]]
    bin.model.acid <- models[[4]]
    bin.model.NN <- models[[5]]
    bin.model.het <- models [[6]]
    bin.model.CO <- models[[8]]
    bin.model.bs <- models[[9]]
    lasso.md.aliph <- models[[7]]
    pls.md.phenolics <- models[[3]]
  }
  mod <- sapply(models, class)
  
  SS <- predict(models[[1]], pX, type="response")
  #acid <- predict(bin.model.acid, pX, type="response")
  NN <- predict(models[[5]], pX, type="response")
  ppX <- as.matrix(pX)
  #dai <- ppls::X2s(ppX, reduce.knots = TRUE)
  #CO <- ppls::new.penalized.pls(b_pPLS.md.CO, dai$Z)
  #CO <- CO$ypred
  CO <- predict(models[[8]], as.matrix(pX))
  bs <- predict(models[[9]], pX, type="response")
  phenolic <- ppls::new.penalized.pls(models[[3]], as.matrix(pX))
  phenolics <- phenolic$ypred
  names(phenolics) <- phenolics
  #aliph <- ppls::new.penalized.pls(b_pPLS.md.aliph, dai$Z)
  #aliph <- aliph$ypred
  aliph <- predict(models[[7]], as.matrix(pX))
  strtn <- cbind(SS, phenolics, NN, aliph, CO, bs)
  names(strtn) <- c("SS", "phenolics", "NN", "aliph", "CO", "bs")
  strtn <- as.data.frame(strtn)
  strtn <- round(strtn)
  strtn[strtn<0] <- 0
  names(strtn)[1] <- "Sulfur"
  names(strtn)[2] <- "phenolics"
  names(strtn)[3] <- "Nitrogen (NN)"
  names(strtn)[4] <- "aliph"
  names(strtn)[5] <- "CO ratio"
  names(strtn)[6] <- "Glycosides (bs)"
  return(strtn)
}

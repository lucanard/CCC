#' @title dataset construction for the application of the CCC method
#' @description the function builds the dataset to apply the the CCC_method to an in-house dataset of standards. The new-built must be saved as model 
#' @param x must be a data.frame containing the compounds' names ("compound"), the monoisotopic mass ("MM"), the retention time ("RT"), the chemical formula ("Formula") and the smiles ("smiles")
#' @usage dataset.building(x)
#' @return the basic dataset to build a CCC_method
#' @export "dataset.building"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples 
#' STD_RP <- read.csv(system.file("extdata", "STD_RP.csv", package = "CCC"), row.names = 1, stringsAsFactors = F)
#' X1Y <- dataset.building(STD_RP)
dataset.building <- function(x) {
  STD_RP <- read.csv(system.file("extdata", "STD_RP.csv", package = "CCC"), row.names = 1, stringsAsFactors = F)
  if (identical(colnames(STD_RP), colnames(x)) == FALSE) {stop("please check the colnames of your data.matrix. The colnames must correspond to the ones in the STD_RP data.frame furnished within the package")}
  options(warnings = F)
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
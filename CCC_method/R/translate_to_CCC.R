CCC_code <- function(x) {
  options(warn = -1)
  mat <- function(x){
    phe <- character(length=length(x))
    for (i in 1:length(x))
      if (x[[i]][1][[1]] == TRUE)
      {phe[i] <- length(x[[i]][2][[1]])} else {phe[i] <- 0}
    return(phe)}
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
  A <- ChemmineR::smiles2sdf(x)
  Pyri = "c1ccncc1"
  Pyrr = "c1cc[nH]c1"
  fur = "c1ccoc1"
  ncyc = "c1cncn1"
  query <- 'c1ccccc1'
  query1 <- "OCC(O)CO"
  query2 <- "O=CO"
  sdfset <- ChemmineR::smiles2sdf(c(Pyri, Pyrr, fur, ncyc))
  sdfset@ID <- c("Pyri", "Pyrr", "fur", "ncyc")
  mcs <- fmcsR::fmcsBatch(sdfset, A, al = 0, au = 0, bl= 0, bu = 0, matching.mode = "aromatic", numParallel = as.numeric(parallel::detectCores()))
  mcs[mcs[,5] < 1] <- 0
  het <- drop(mcs[,5])
  mols <- tryCatch({sapply(x, rcdk::parse.smiles)}, error = function (e) {sapply(x, function(x) rcdk::parse.smiles(x, kekulise = FALSE))})
  #does2 <- rcdk::matches(query2, mols, return.matches = FALSE)
  does <- rcdk::matches(query, mols, return.matches = TRUE)
  does1 <- rcdk::matches(query1, mols, return.matches = FALSE)
  #acidic <- as.numeric(does2 + 0)
  bs <- as.numeric(does1 + 0)
  phenolics <- as.numeric(mat(does))
  Y <- cbind(sulfur, phenolics, NN, aliph, CO, bs)
  return(Y)
}
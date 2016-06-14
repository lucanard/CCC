autoranking.metfRag <- function(score, sorts, queries, i2) {
  if (length(score) == 0) {return(score)} else {
  files <- score
  max.peaks <- as.numeric(max(files$NoExplPeaks))
  files <- files[as.numeric(files$NoExplPeaks) >= max.peaks/10,]
  com <- CCC_code(as.character(files$SMILES))
  predicted <- sorts[sorts$rowname == queries$name[i2], 10:16]
  ebb <- apply(com, 1, '-', as.numeric(predicted))
  eb <- colSums(abs(ebb))
  eb <- eb + 1
  reordi <- order(eb)
  candidate <- files[reordi,]
  print("candidates ranked according to CCC method")
  }
  return(candidate)
}
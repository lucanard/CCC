off.ranking.metfRag <- function(scores, tni, marker) {
  compi <- do.call(rbind, tni)
  files <- scores[[1]]
  max.peaks <- as.numeric(max(files$NoExplPeaks))
  files <- files[as.numeric(files$NoExplPeaks) >= max.peaks/10,]
  com <- CCC_code(as.character(files$SMILES))
  predicted <- compi[compi$rowname == names(score), 10:16]
  ebb <- apply(com, 1, '-', as.numeric(predicted))
  eb <- colSums(abs(ebb))
  eb <- eb + 1
  reordi <- order(eb)
  candidate <- files[reordi,]
  names(candidate) <- names(score)
  return(candidate)
}


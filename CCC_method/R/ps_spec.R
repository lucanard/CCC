#' @title Pseudospectra extraction for MSe data analysis with metfRag
#' @description The function extracts all the pseudospectra given by the CAMERA grouping, which also have at least one isotopic group in the grouped isotopic pattern given by the apply.model function 
#' @param "tni: the list of grouped features obtained from the function apply.model" 
#' @param "peakTable: the peakTable obtained from the runLC function from metaMS" 
#' @usage ps_spec(tni, peakTable)
#' @return "ps_spectra: all the pseudospectra per each group of feature from the list tni"
#' @export "ps_spec"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples 
#' data("PeakTable")
#' tni <- apply.model(peakTable, polarity = "negative")
#' ps_spectra <- ps_spec(tni, peakTable)
ps_spec <- function(tni, peakTable) {
  compi <- do.call(rbind, tni)
  groups <- unique(compi$gro)
  mzs <- unique(peakTable$mz[peakTable$pcgroup %in% groups])
  int <- rowSums(peakTable[peakTable$pcgroup %in% groups, 8:ncol(peakTable)])
  gro <- peakTable$pcgroup[peakTable$pcgroup %in% groups]
  ps <- data.frame(mzs, int, gro)
  return(ps)
}

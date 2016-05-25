#' @title Pseudospectra extraction for MSe data analysis with metfRag
#' @description The function extracts all the pseudospectra given by the CAMERA grouping, which also have at least one isotopic group in the grouped isotopic pattern given by the apply.model function 
#' @param tni the list of grouped features obtained from the function apply.model
#' @param peakTable the peakTable obtained from the runLC function from metaMS
#' @param peaklist a peaklist object obtain from the function getpeaklist from CAMERA 
#' @usage ps_spec(tni, peakTable, peaklist)
#' @return "ps_spectra: all the pseudospectra per each group of feature from the list tni"
#' @export "ps_spec"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples 
#' data("PeakTable")
#' tni <- apply.model(peakTable, polarity = "negative")
#' ps_spectra <- ps_spec(tni, peakTable)
ps_spec <- function(tni, peakTable = NULL, peaklist = NULL) {
  compi <- do.call(rbind, tni)
  groups <- unique(compi$gro)
  if (!is.null(peakTable)  & !is.null(peaklist)) {stop("load only a peakTable or a peaklist, avoid to load both")}
  if (is.null(peakTable) & is.null(peaklist)) {stop("load at least a peakTable or a peaklist")}
  if (!is.null(peakTable)) {
  mzs <- unique(peakTable$mz[peakTable$pcgroup %in% groups])
  int <- rowSums(peakTable[peakTable$pcgroup %in% groups, 8:ncol(peakTable)])
  gro <- peakTable$pcgroup[peakTable$pcgroup %in% groups]
  ps <- data.frame(mzs, int, gro)
  }
  if (!is.null(peaklist)) {
    mzs <- unique(peaklist$mz[peaklist[,ncol(peaklist)] %in% groups])
    int <- rowSums(peaklist[peaklist[,ncol(peaklist)] %in% groups, 9:ncol(peaklist)])
    gro <- peaklist[,ncol(peaklist)][peaklist[,ncol(peaklist)] %in% groups]
    ps <- data.frame(mzs, int, gro)
  }
  return(ps)
}

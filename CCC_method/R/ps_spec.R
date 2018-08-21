#' @title Pseudospectra extraction for MSe data analysis with metfRag
#' @description The function extracts all the pseudospectra given by the CAMERA grouping, which also have at least one isotopic group in the grouped isotopic pattern given by the apply.model function 
#' @param tni the list of grouped features obtained from the function apply.model
#' @param peakTable the peakTable obtained from the runLC function from metaMS
#' @param xsaFA the xsaFA object obtained from CAMERA 
#' @usage ps_spec(tni, peakTable = NULL, xsaFA = NULL)
#' @return ps_spectra all the pseudospectra per each group of feature from the list tni
#' @import CAMERA
#' @export "ps_spec"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples 
#' data("peakTable")
#' tni <- CCC.peakTable(peakTable, polarity = "negative")
#' ps_spectra <- ps_spec(tni, peakTable)
ps_spec <- function(tni, peakTable = NULL, xsaFA = NULL) {
  compi <- do.call(rbind, tni)
  groups <- unique(compi$gro)
  if (!is.null(peakTable)  & !is.null(xsaFA)) {stop("load only a peakTable or a xsaFA object, avoid to load both")}
  if (is.null(peakTable) & is.null(xsaFA)) {stop("load at least a peakTable or a xsaFA object")}
  if (!is.null(peakTable)) {
  mzs <- unique(peakTable$mz[peakTable$pcgroup %in% groups])
  int <- rowSums(peakTable[peakTable$pcgroup %in% groups, 8:ncol(peakTable)])
  gro <- as.numeric(peakTable$pcgroup[peakTable$pcgroup %in% groups])
  ps <- data.frame(mzs, int, gro)
  ps <- ps[order(ps$gro),]
  }
  if (!is.null(xsaFA)) {
    peaklist <- CAMERA::getPeaklist(xsaFA)
    mzs <- unique(peaklist$mz[peaklist$pcgroup %in% groups])
    int <- rowSums(peaklist[peaklist$pcgroup %in% groups, 9:(ncol(peaklist)-4)])
    gro <- as.numeric(peaklist$pcgroup[peaklist$pcgroup %in% groups])
    ps <- data.frame(mzs, int, gro)
    ps <- ps[order(ps$gro),]
  }
  return(ps)
}

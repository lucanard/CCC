#' @title Pseudospectra extraction for MSe data analysis with metfRag
#' @description The function extracts all the pseudospectra given by the CAMERA grouping, which also have at least one isotopic group in the grouped isotopic pattern given by the apply.model function 
#' @param "tni: the list of grouped features obtained from the function apply.model" 
#' @param "peakTable: the peakTable obtained from the runLC function from metaMS" 
#' @usage ps_spec(tni, peakTable)
#' @return "ps_spectra: all the pseudospectra per each group of feature from the list tni"
#' @export "ps_spec"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples 
#' tni <- apply.model(peakTable, polarity = "negative")
#' ps_spectra <- ps_spec(tni, peakTable)
ps_spec <- function(tni, peakTable){
  compi <- do.call(rbind, tni)
  row.names(compi) <- compi$rowname
  psspectra <- list()
  for (i in 1:length(compi$gro)) {
    mzo = peakTable$mz[which(peakTable$pcgroup == compi$gro [i])]
    int = rowSums(peakTable[which(peakTable$pcgroup == compi$gro [i]),8:(ncol(peakTable))])
    temp1 = data.frame(mzo, int)
    psspectra[[i]] = drop(temp1)
  }
  names(psspectra) <- compi$gro
  order <- order(names(psspectra))
  psspectra <- unique(psspectra[order])
  return(psspectra)
}

#' @title Ranking function to get the best candidates from the metfRag analysis from .sdf files
#' @description The function reorders the metfrag candidates using the information given by the CCC_method. It uploads the .sdf files from a directory, and performs reordering of the data inside each file. It returns a list of SDF files each containing a list of candidates. NB: the .sdf file names must correspond to the original raw number of the markers in the pekaTable. 
#' @param directory where the sdf files are stored. The file name must correspond to the feature number (rowname) in the original peaktable 
#' @param tni a list of grouped features given by the apply.model function 
#' @usage ranking.2016(directory, tni)
#' @return a list of SDF files reordered according to the CCC method parameters
#' @export "ranking.2016"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples
#' data (peakTable)
#' tni <- CCC.peakTable(peakTable, polarity = "negative")
#' candidates <- ranking.2016(system.file("extdata", "2016", package = "CCC"), tni)
ranking.2016 <- function(directory, tni) {
  setwd(directory)
  compi <- do.call(rbind, tni)
  row.names(compi) <- compi$rowname
  file <- list.files(path = ".", full.names = FALSE)
  fi <- str_detect(file, pattern = "csv")  
  fil <- str_replace(file[fi], pattern = ".csv", replacement="")
  fil <- as.numeric(fil)
  candidate <- list()
  for (i in 1:length(fil)) {
    print(fil[i])
    files <- read.csv(file[fi][i], stringsAsFactors=FALSE)
    max.peaks <- as.numeric(max(files$NoExplPeaks))
    files <- files[as.numeric(files$NoExplPeaks) >= max.peaks/10,]
    es <- parse.inchi(files$InChI)
    esi <- sapply(es, function(x) length(x))
    es <- es[esi != 0]
    esi <- sapply(es, function(x) get.smiles(x))
    com <- CCC_code(esi)
    predicted <- compi[compi$rowname == fil[i], 10:16]
    ebb <- apply(com, 1, '-', as.numeric(predicted))
    eb <- colSums(abs(ebb))
    eb <- eb + 1
    reordi <- order(eb)
    correct <- files[reordi,]
    candidate[[i]] <- correct
    names(candidate[[i]]) <- names(files)
  }
  return(candidate)
}
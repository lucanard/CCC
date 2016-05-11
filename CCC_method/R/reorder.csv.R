#' Reordering function to get the best candidates from the metfRag analysis from .csv files
#'
#' @param directory where the .csv files are stored. The file name must correspond to the feature number (rowname) in the original peaktable  
#' @param tni: a list of grouped features given by the apply.model function 
#'
#' @return a list of csv files reordered according to the CCC method parameters
#' @export "reorder.csv"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples
reorder.csv <- function(directory, tni) {
  setwd(directory)
  compi <- do.call(rbind, tni)
  row.names(compi) <- compi$rowname
  file <- list.files(pattern = "csv", full.names = FALSE)
  fil <- str_replace(file, pattern = ".csv", replacement="")
  fil <- as.numeric(fil)
  candidate <- list()
  for (i in 1:length(fil)) {
    print(fil[i])
    files <- read.csv(file[i], stringsAsFactors=FALSE)
    name <- files[4,]
    files <- files[5:nrow(files),]
    names(files) <- name
    es <- files$Smiles
    com <- CCC_code(es)
    predicted <- compi[row.names(compi) == fil[i], 9:15]
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
#' Reordering function to get the best candidates from the metfRag analysis from .csv files
#'
#' @param directory where the .csv files are stored 
#' @param compi: a data.frame containing all the features extracted by the apply.model function 
#'
#' @return a list of csv files reordered according to the CCC method parameters
#' @export "reorder.csv"
#'
#' @examples
reorder.csv <- function(directory, compi) {
  setwd(directory)
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
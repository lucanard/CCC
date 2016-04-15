reorder.csv <- function(directory, compi, selection) {
  setwd(directory)
  file <- list.files(pattern = "csv", full.names = FALSE)
  candidate <- list()
  for (i in 1:length(file)) {
    print(file[i])
    files <- read.csv(file[i], stringsAsFactors=FALSE)
    name <- files[4,]
    files <- files[5:nrow(files),]
    names(files) <- name
    es <- files$Smiles
    com <- CCC_code(es)
    predicted <- compi[selection, 10:16]
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
reorder.sdf <- function(directory, compi, selection) {
  setwd(directory)
  file <- list.files(pattern = "sdf", full.names = FALSE)
  fil <- str_split(file, pattern = "_")
  fi <- do.call(cbind, fil)
  fil <- fi[1,]
  candidate <- list()
  for (i in 1:length(file)) {
    print(fil[i])
    esp <- read.SDFstr(file[i])
    esp <- as(esp, "SDFset")
    espo <- sdf2smiles(esp)
    es <- toString(espo)
    es <- as.vector(unlist(strsplit(es, split=", ", perl = TRUE)))
    com <- CCC_code(es)
    predicted <- compi[selection, 10:16]
    ebb <- apply(com, 1, '-', as.numeric(predicted))
    eb <- colSums(abs(ebb))
    eb <- eb + 1
    reordi <- order(eb)
    correct <- esp[reordi,]
    candidate[[i]] <- correct
  }
  return(candidate)
}
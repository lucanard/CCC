#' Reordering function to get the best candidates from the metfRag analysis from .sdf files
#'
#' @param directory where the sdf files are stored 
#' @param compi: a data.frame containing all the features extracted by the apply.model function
#'
#' @return a list of SDF files reordered according to the CCC method parameters
#' @export "reorder.sdf"
#'
#' @examples 
reorder.sdf <- function(directory, compi) {
  setwd(directory)
  file <- list.files(pattern = "sdf", full.names = FALSE)
  fil <- str_replace(file, pattern = ".sdf", replacement="")
  fil <- as.numeric(fil)
  candidate <- list()
  for (i in 1:length(fil)) {
    print(fil[i])
    esp <- read.SDFstr(file[i])
    esp <- as(esp, "SDFset")
    espo <- sdf2smiles(esp)
    es <- toString(espo)
    es <- as.vector(unlist(strsplit(es, split=", ", perl = TRUE)))
    com <- CCC_code(es)
    predicted <- compi[row.names(compi) == fil[i], 9:15]
    ebb <- apply(com, 1, '-', as.numeric(predicted))
    eb <- colSums(abs(ebb))
    eb <- eb + 1
    reordi <- order(eb)
    correct <- esp[reordi,]
    candidate[[i]] <- correct
  }
  return(candidate)
}
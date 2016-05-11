#' @title Reordering function to get the best candidates from the metfRag analysis from .sdf files
#' @description 
#' @param directory where the sdf files are stored. The file name must correspond to the feature number (rowname) in the original peaktable 
#' @param tni: a list of grouped features given by the apply.model function 
#'
#' @return a list of SDF files reordered according to the CCC method parameters
#' @export "reorder.sdf"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples 
reorder.sdf <- function(directory, tni) {
  setwd(directory)
  compi <- do.call(rbind, tni)
  row.names(compi) <- compi$rowname
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
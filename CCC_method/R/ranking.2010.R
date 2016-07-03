#' @title Ranking function to get the best candidates from the metfRag analysis from .csv files
#' @description The function reorders the metfrag candidates using the information given by the CCC_method. It uploads the .csv files from a directory, and performs reordering of the data inside each file. It returns a data.frame containing all the files together. NB: the .csv files must be from the web metfrag2010 and must be already converted in csv files (not xls). The name of the .csv file must correspond to the original raw number of the markers in the peakTable. 
#' @param directory where the .csv files are stored. The file name must correspond to the feature number (rowname) in the original peaktable  
#' @param tni is a list of grouped features given by the CCC.peakTable/CCC.xsaFA function
#' @param sep is the separator preferred in the kind of file (useful only for csv files). Default is the "comma". Any separator can be applied. 
#' @usage ranking.2010(directory, tni, sep)
#' @return a list of csv files reordered according to the CCC method parameters
#' @export "ranking.2010"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples
#' data (peakTable)
#' tni <- CCC.peakTable(peakTable, polarity = "negative")
#' candidates <- ranking.2010(system.file("extdata", "2010", package = "CCC"), tni)
ranking.2010 <- function(directory, tni, sep) {
  setwd(directory)
  if (all(str_detect(file, pattern = "csv") != TRUE) & all(str_detect(file, pattern = "sdf") != TRUE)) {
    stop("csv/sdf files are missing: please control the directory")
  }
  if (any(str_detect(file, pattern = "csv") == TRUE) & any(str_detect(file, pattern = "sdf") == TRUE)) {
    warning("both csv and sdf files are present in the folder only the former will be used")
  }
  compi <- do.call(rbind, tni)
  row.names(compi) <- compi$rowname
  file <- list.files(path = ".", full.names = FALSE)
  if (any(str_detect(file, pattern = "csv") == TRUE)) {
  fi <- str_detect(file, pattern = "csv")  
  fil <- str_replace(file[fi], pattern = ".csv", replacement="")
  fil <- as.numeric(fil)
  candidate <- list()
  for (i in 1:length(fil)) {
    print(fil[i])
    files <- read.csv(file[fi][i], sep = ",", stringsAsFactors=FALSE)
    if (length(files) == 1) {stop("please select the prope separator for your csv files")}
    name <- files[4,]
    files <- files[5:nrow(files),]
    names(files) <- name
    max.peaks <- as.numeric(max(files$`# of Peaks Explained`))
    files <- files[as.numeric(files$`# of Peaks Explained`) >= max.peaks/10,]
    es <- as.character(files[["Smiles"]])
    com <- CCC_code(es)
    predicted <- compi[compi$rowname == fil[i], 10:16]
    ebb <- apply(com, 1, '-', as.numeric(predicted))
    eb <- colSums(abs(ebb))
    eb <- eb + 1
    reordi <- order(eb)
    correct <- files[reordi,]
    candidate[[i]] <- list(correct, files)
    }
  } else {
    fo <- str_detect(file, pattern = "sdf")
    fil <- str_replace(file[fo], pattern = ".sdf", replacement="")
    fil <- as.numeric(fil)
    candidate <- list()
    for (i in 1:length(fil)) {
      print(fil[i])
      esp <- read.SDFstr(file[fo][i])
      esp <- as(esp, "SDFset")
      selic <- numeric()
      for (i2 in 1:length(esp)) { 
        selic[i2] <- as.numeric(drop(esp[[i2]][[4]][[5]]))
      }
      esp <- esp[selic >= max(selic/10)]
      espo <- sdf2smiles(esp)
      es <- toString(espo)
      es <- as.vector(unlist(strsplit(es, split=", ", perl = TRUE)))
      com <- CCC_code(es)
      predicted <- compi[compi$rowname == fil[i], 10:16]
      ebb <- apply(com, 1, '-', as.numeric(predicted))
      eb <- colSums(abs(ebb))
      eb <- eb + 1
      reordi <- order(eb)
      correct <- esp[reordi,]
      candidate[[i]] <- list(correct, esp)
    }
  }
  return(candidate)
}
#' @title Wrapper function to online querying of m/z candidates in metfRag and get CCC re-ranking
#' @description The function handles MSe data to obtain MS/MS analysis through the metfRag package and to reorder the molecules according to the CCC method calculations. It directly queries the the 4 most intense ions of a pseudospectra in metfRag, discards imporbabloe structures according to some heuristic rules, then reranks the obtained candidates using the Y variables extrapolated in the reduced peakTable tni.
#' @param tni the reduced peakTable given by the CCC.peaktable or CCC.xsaFA functions
#' @param ps_spectra the pseudospectra given by the ps_spec function
#' @param threshold the threshold used to discard low intense MS/MS markers
#' @param markers a list of markers given from your own statistical analysis. If NULL, all the first most intense ion of each speudospectra are analyzed; WARNING = data file must be massive and it might lead to errors.
#' @param pre_filter if TRUE, it uses the features measured by the CCC method to pre-filter candidates before performing in-silico fragmentation. It speeds up the process, but it is error-prone. default is set to FALSE
#' @param database the database to use to fetch candidate structures with no default. Please choose between "PubChem", "ExtendedPubChem", "ChemSpider" (it requires Chemspider token), "KEGG".
#' @param token your personal ChemSpider token
#' @param n_candidates the number of candidates parent ions to be considered fro MS/MS analysis from the same spectrum. 
#' @usage metfragging(tni, ps_spectra, threshold, markers, pre_filter, n_candidates, database, token)
#' @return a list of data.frames containing the outcome of metfRag reordered according to the CCC method
#' @export "metfragging"
#' @author Luca Narduzzi "nardluca@gmail.com"
#' @examples 
#' data(peakTable)
#' tni <- CCC.peakTable(peakTable, polarity = "negative")
#' ps_spectra <- ps_spec(tni, peakTable)
#' results <- metfragging(tni, ps_spectra, markers = 2713, database = "KEGG")
metfragging <- function(tni, ps_spectra, threshold = 1000, markers = NULL, pre_filter = FALSE, n_candidates = 2, database, token = NULL) {
  if (missing(tni) == TRUE) {stop("please insert reduced peakTable")}
  if (missing(ps_spectra) == TRUE) {stop("please insert pseudo-spectras")}
  if (missing(database) == TRUE) {stop("please insert the database to fetch the structures")}
  if (pre_filter == TRUE) {warning("the pre-filter function is still not complete and it might lead to errors")}
  if (database == "ChemSpider" & is.null(token)) {stop("please insert ChemSpider token")}
  if (database != "ChemSpider" & !is.null(token)) {token <- NULL}
  score <- list()
scores <- list()
compi <- do.call(rbind, tni)
if (!is.null(markers)) {groups <- unique(compi$gro[compi$rowname %in% markers])} else {groups <- names(tni)}
for (i in 1:length(groups)) {
  sorti <- compi[compi$gro == groups[i],]
    sorts <- sorti[order(sapply(sorti$int, max), decreasing=TRUE),]
    mark <- which(sorts$rowname %in% markers)
    name <- as.numeric(sorts$rowname)
    thr <- as.numeric(sorts$mzs)
    if (length(thr) >= n_candidates) {thr <- thr[1:n_candidates]} else {thr <- thr[1:max(length(thr))]}
    if (length(name) >= n_candidates) {name <- name[1:n_candidates]} else {name <- name[1:max(length(name))]}
    queries <- data.frame(name, thr)
    if (any(as.logical(name)) %in% sorts[mark,"rowname"]) {queries <- queries} else {thr1 <- sorts[mark, "mzs"]
    name1 <- as.character(sorts[mark,"rowname"])
    queries1 <- data.frame(name1, thr1)
    names(queries1) <- c("name", "thr")
    queries <- unique(rbind(queries, queries1))
    }
    for (i2 in 1:nrow(queries)) {
    peaklist <- as.matrix(ps_spectra[which(ps_spectra$gro == groups[i]), -3])
    if (pre_filter == TRUE) {SMARTS <- code_CCC(sorts[i2,10:16])
      filters <- "SmartsSubstructureInclusionScore"
      met_score <- c(1.0, 1.0)} else {SMARTS <- NULL
        filters <- NULL
        met_score <- NULL}
    maxnC <- as.numeric(round(queries$thr[i2]/12))
    minnC <- as.numeric(round(maxnC/8))
    maxnH <- as.numeric(round(maxnC*2))
    minnH <- minnC
    if (minnH <= 1) {minnH <- 2} else {minnH <- minnH}
    if (sorts$CO[i2] == 1 & pre_filter == TRUE) {maxnO = round(maxnC/7)} else {maxnO <- round(maxnC * 1.2)}
    minnO <- round(minnC*0.2)
    maxnN <- 10
    minnN <- 0
    maxnP <- 3
    minnP <- 0
    maxnS <- 3
    minnS <- 0
    limitup <- str_c(c("C", maxnC, "H", maxnH, "N", maxnN, "O",maxnO, "P", maxnP, "S", maxnS), collapse = "")
    limitdown <- str_c(c("C", minnC, "H", minnH, "N",minnN, "O", minnO, "P", minnP, "S", minnS), collapse = "")
    settingsObject<-list()
    settingsObject[["DatabaseSearchRelativeMassDeviation"]] <- 10.0
    settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]] <- 0.01
    settingsObject[["FragmentPeakMatchRelativeMassDeviation"]]<-10.0
    settingsObject[["MinimumAbsolutePeakIntensity"]] <- threshold
    settingsObject[["NumberThreads"]] <- as.numeric(parallel::detectCores())
    settingsObject[["MetFragDatabaseType"]] <- database
    settingsObject[["MetFragPreProcessingCandidateFilter"]] <- c("UnconnectedCompoundFilter","IsotopeFilter", "MinimumElementsFilter", "MaximumElementsFilter", "ElementExclusionFilter")  
    settingsObject[["MetFragPostProcessingCandidateFilter"]]<- "InChIKeyFilter"
    settingsObject[["NeutralPrecursorMass"]]<- queries$thr[i2]
    settingsObject[["PeakList"]] <- peaklist
    settingsObject[["ScoreSmartsInclusionList"]]<- c(unlist(SMARTS))
    settingsObject[["FilterMinimumElements"]] <- limitdown
    settingsObject[["FilterMaximumElements"]] <- limitup
    settingsObject[["ChemSpiderToken"]] <- token
    settingsObject[["FilterExcludedElements"]] <- c("Si", "Br", "Fe", "Mn", "K", "F", "Na", "Cl", "Mg", "Ca", "Se", "I", "Co", "Ni", "Cr", "Be")
    #settingsObject[["PrecursorIonMode"]] <- 0
    settingsObject[["MetFragScoreTypes"]] <- c("FragmenterScore", filters)
    settingsObject[["MetFragScoreWeights"]]<- met_score
    candidate <- metfRag::run.metfrag(settingsObject)
    if (pre_filter == FALSE) {score[[i2]] <- autoranking.metfRag(as.data.frame(candidate), sorts, queries, i2)} else {score[[i2]] = candidate}
  }
    names(score) <- name
    scores[[i]] <- score
  }
  return(scores)
}

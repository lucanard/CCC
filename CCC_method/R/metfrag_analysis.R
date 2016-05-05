metfragging <- function(tni, ps_spectra, markers) {
score <- list()
scores <- list()
if (pol == "positive") {mode <- 1} else {mode <- -1}
compi <- do.call(rbind, tni)
compx <- compi[compi$rowname %in% markers,]
groups <- compx$gro
for (i in groups) {
  if (length(tni[[i]]$mz)  >= 2) {
    sort <- tni[[i]][order(sapply(tni[[i]]$int, function(x) x[1], simplify=TRUE), decreasing=TRUE),]
    thr <- as.numeric(sort$mzs)} else
    {thr <- tni[[i]]$mzs[1]}
  if (length(thr) >= 5) {thr <- thr[1:4]} else {thr <- thr[1:max(length(thr))]}
  for (i2 in 1:length(thr)) {
    settingsObject<-list()
    settingsObject[["DatabaseSearchRelativeMassDeviation"]]<-10.0
    settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]]<-0.005
    settingsObject[["FragmentPeakMatchRelativeMassDeviation"]]<-10.0
    #settingsObject[["MinimumAbsolutePeakIntensity"]] <- as.numeric(quantile(ps_spectra[[i]][,2], probs = 0.25)[[1]])
    #settingsObject[["NumberThreads"]] <- 2
    settingsObject[["IsPositiveIonMode"]] <- FALSE
    settingsObject[["MetFragDatabaseType"]] <- "PubChem"
    #settingsObject[["MetFragPreProcessingCandidateFilter"]] <- c("ElementExclusionFilter", "ElementInclusionFilter", "IsotopeFilter", "MaximumElementsFilter", "MinimumElementsFilter", "UnconnectedCompoundFilter")  
    #settingsObject[["ChemSpiderToken"]] <- "8a8a28e1-a41e-489c-907b-ac8c8a6f1a76"
    settingsObject[["NeutralPrecursorMass"]]<- thr[i2]
    #settingsObject[["MaximumTreeDepth"]] <- c(1, 2)
    settingsObject[["PeakList"]] <- as.matrix(ps_spectra[[i]])
    settingsObject[["FilterMinimumElements"]] <- "C2H4N0O1P0S0"
    settingsObject[["FilterMaximumElements"]] <- "C50H100N10O50P3S3"
    settingsObject[["FilterIncludedElements"]] <- list("C", "H", "N", "O", "P", "S")
    settingsObject[["FilterExcludedElements"]] <- list("Si", "Br", "Fe", "Mn", "K", "F")
    #settingsObject[["PrecursorIonMode"]] <- 0
    #settingsObject[["MetFragScoreTypes"]] <- c("FragmenterScore", "RetentionTimeScore", "CombinedReferenceScore")
    score[[i2]] <- run.metfrag(settingsObject)
    }
  scores[[i]] <- score[[i2]]
  }
  return(scores)
}
metfragging <- function(tni, ps_spectra, threshold = 1000, markers = NULL) {
score <- list()
scores <- list()
compi <- do.call(rbind, tni)
if (!is.null(markers)) {groups <- unique(compi$gro[compi$rowname == markers])} else {groups <- names(tni)}
for (i in 1:length(groups)) {
  sorti <- compi[compi$gro == groups[i],]
    sort <- sorti[order(sapply(sorti$int, max), decreasing=TRUE),]
    thr <- as.numeric(sort$mzs)
    if (length(thr) >= 5) {thr <- thr[1:4]} else {thr <- thr[1:max(length(thr))]}
  for (i2 in 1:length(thr)) {
    peaklist <- as.matrix(ps_spectra[which(ps_spectra$gro == groups[i]), -3])
    settingsObject<-list()
    settingsObject[["DatabaseSearchRelativeMassDeviation"]]<-10.0
    settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]]<-0.005
    settingsObject[["FragmentPeakMatchRelativeMassDeviation"]]<-10.0
    settingsObject[["MinimumAbsolutePeakIntensity"]] <- threshold
    settingsObject[["NumberThreads"]] <- 2
    settingsObject[["MetFragDatabaseType"]] <- "PubChem"
    settingsObject[["MetFragPreProcessingCandidateFilter"]] <- c("UnconnectedCompoundFilter","IsotopeFilter", "MinimumElementsFilter", "MaximumElementsFilter", "ElementInclusionFilter", "ElementExclusionFilter")  
    settingsObject[["MetFragPostProcessingCandidateFilter"]]<- "InChIKeyFilter"
    #settingsObject[["ChemSpiderToken"]] <- "8a8a28e1-a41e-489c-907b-ac8c8a6f1a76"
    settingsObject[["NeutralPrecursorMass"]]<- thr[i2]
    #settingsObject[["MaximumTreeDepth"]] <- c(1, 2)
    settingsObject[["PeakList"]] <- peaklist
    settingsObject[["FilterMinimumElements"]] <- "C2H4N0O1P0S0"
    settingsObject[["FilterMaximumElements"]] <- "C50H100N10O50P3S3"
    settingsObject[["FilterIncludedElements"]] <- c("C", "H", "O", "N", "S", "P")
    settingsObject[["FilterExcludedElements"]] <- c("Si", "Br", "Fe", "Mn", "K", "F", "Na", "Cl")
    #settingsObject[["PrecursorIonMode"]] <- 0
    #settingsObject[["MetFragScoreTypes"]] <- c("FragmenterScore", "RetentionTimeScore", "CombinedReferenceScore")
    score[[i2]] <- run.metfrag(settingsObject)
    }
  scores[[i]] <- score[[i2]]
  }
  return(scores)
}




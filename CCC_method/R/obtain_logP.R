require(rcdk)
logping <- function(x)
molecules <- tryCatch({sapply(x, rcdk::parse.smiles)}, error = function (e) {sapply(x, function(x) rcdk::parse.smiles(x, kekulise = FALSE))})
xlogp <- sapply(molecules, function(x) get.xlogp(x))
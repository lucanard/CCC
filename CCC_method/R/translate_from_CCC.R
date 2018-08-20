code_CCC <- function(x) {
 smarties <- list()
 if (x$SS == 1) {smarties[[1]] <- "S"} else {smarties[[1]] <- NULL}
 if (x$phenolics >= 1) {smarties[[2]] <- rep("c1ccccc1", x$phenolics)} else {smarties[[2]] <- NULL}
 #if (x$acid == 1) {smarties[[3]] <- "O=CO"} else {smarties[[3]] <- NULL}
 if (x$NN == 1)  {smarties[[4]] <- "N"} else {smarties[[3]] <- NULL}
 if (x$aliph >= 1) {smarties[[5]] <- "CCCC"} else {smarties[[4]] <- NULL}
 if (x$bs == 1) {smarties[[5]] <- "OCC(O)CO"} else {smarties[[6]] <- NULL}
return(smarties)
}
  
  
  
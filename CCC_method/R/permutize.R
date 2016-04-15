permutize <- function (x, per.test = FALSE){
  if (per.test != FALSE){
    x <- lapply(x, function(x) x[sample(length(x))])
  } else {x <- x}
  return(x)
}
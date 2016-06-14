permutize <- function (x){
  x <- lapply(x, function(x) x[sample(length(x))])
  return(x)
}
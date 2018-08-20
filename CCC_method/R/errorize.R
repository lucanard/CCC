errorize <- function (x, errori, i2) {
  ppms <- c(-50, -30, -10, -5, -3, 0, 3, 5, 10, 30, 50)
  Cees <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5)
  X2 <- data.frame()
  Xm <- x
  if (errori == "ppm") {
    x[,2] <- (x[,2]+(x[,2]/1000000*ppms[i2]))
    x[,4] <- (x[,2] - round(x[,2]))
    x[,5] <- ((x[,4]/x[,2])*1000000)
    x[,7] <- ((x[,2]-(x[,3]*12) - round(x[,2]-(x[,3]*12)))/round(x[,2]-(x[,3]*12)))*1000000
  }
  if (errori == "Cees") {
    x[,3] <- (x[,3] + Cees[i2])
    maxC <- round(x[,2]/12)
    wC <- which(x[,3] - maxC >= 0)
    x[wC,3] <- maxC[wC] - 1
    x[which(x[,3] <= 0), 3] <- 1
    x[,6] <- (x[,3]*12/x[,2])*100
    x[,7] <- ((x[,2]-(x[,3]*12) - round(x[,2]-(x[,3]*12)))/round(x[,2]-(x[,3]*12)))*1000000
  }
  if (errori == "nerr") {
    x[,2] <- (x[,2]+rnorm(x[,2], mean=0, sd=0.001))
    x[,4] <- (x[,2] - round(x[,2]))
    x[,5] <- ((x[,4]/x[,2])*1000000)
    x[,3] <- round(x[,3]+rnorm(x[,3], mean=0, sd=1))
    maxC <- round(x[,2]/12)
    wC <- which(x[,3] - maxC >= 0)
    x[wC,3] <- maxC[wC] - 1
    x[which(x[,3] <= 0), 3] <- 1
    x[,6] <- ((x[,3]*12/x[,2])*100)
    x[,7] <- ((x[,2]-(x[,3]*12) - round(x[,2]-(x[,3]*12)))/(round(x[,2]-(x[,3]*12))))*1000000
  }
  if (errori == "none") {x <- x}
  tab <- data.frame(i2, 11)
  names(tab) <- c("n.test", "tot.test")
  print(tab, row.names = F)
  X2 <- x
  return(X2)
}
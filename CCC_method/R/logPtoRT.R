require(rcdk)
data(FinSTD)
molecules <- parse.smiles(FinSTD$smiles, kekulise=TRUE)
xlogp <- sapply(molecules, function(x) get.xlogp(x))
rts <- FinSTD$RT[which(!is.na(FinSTD$RT))]
xlogp.real <- xlogp[which(!is.na(FinSTD$RT))]
xlogp.pred <- xlogp[which(is.na(FinSTD$RT))]
names(xlogp.real) <- NULL
names(xlogp.pred) <- NULL

plot(xlogp.real, rts)
logpx <- drop(xlogp.real)
rtsy <- rts
fit <- lm(rtsy ~ poly(logpx, 2))

plot(rtsy~logpx)
lines(sort(logpx), fitted(fit)[order(logpx)], col='red', type='b') 
logpx <- data.frame(xlogp.pred)
names(logpx) <- "logpx"
RTs <- predict(fit, newdata = logpx)

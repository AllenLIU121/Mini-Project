library(rlist)
library(plyr)

rouletteSpins <- function(p){
  	x = 1
  	while (runif(1, 0, 1) < p){
      x <- x+1
  	}
  	return (x)
}

p <- 18/37
xGrid <- seq(1,7,1)
N <- 10**5

mc <-vector()
for (i in (1:N)) {
  mc <- append(mc,rouletteSpins(p))
}

mcEstimate <-vector()
for (j in xGrid){
  mcEstimate <- append(mcEstimate,sum(mc==j)/N)
}

gPmf <- dgeom(x = xGrid-1, prob = p)

plot(x = c(xGrid,xGrid), y = c(mcEstimate,gPmf), col = c("blue", "red"), xlab="X", ylab="Probability",pch = c(19,4), cex = 1.5)
lines(x = c(xGrid,xGrid), y = c(mcEstimate,gPmf), col = c("blue", "red"), type = 'h',lwd = 2)
legend(5,0.5,legend=c("mcEstimate", "PMF"),col=c("blue","red"),,pch = c(19,4) , cex=1)




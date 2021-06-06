library(stats)


muRange <- seq(51,55,0.02)
n <- 20
N <- 10**3
mu0 <- 53.0
powerT <- vector()
powerU <- vector()

for (muActual in muRange) {
  rejectT <- 0
  rejectU <- 0
  set.seed(1)
  
  for (i in 1:N) {
    data <- rnorm(n, mean <- muActual,sd <- 1.2)
    xBar <- mean(data)
    stdDev <- sd(data)
    
    tStatT <- (xBar - mu0)/(stdDev/sqrt(n))
    pValT <- 2*(1-pt(abs(tStatT),n-1 ))	
    
    xPositive <- sum(data > mu0)
    uStat <- max(xPositive, n-xPositive)
    pValSign <- 2*(dbinom(uStat,n,0.5))
    rejectT = rejectT + (pValT < 0.05)
    rejectU = rejectU + (pValSign < 0.05)
  }
  
  powerT <- append(powerT,rejectT/N)
  powerU <- append(powerU,rejectU/N)
}

plot(muRange, powerT, type = "l",col="blue", xlim=c(51,55),ylim = c(0, 1),xlab="Different values of muActual",
     ylab="Proportion of times H0 rejected")
lines(muRange, powerU, col="red")
legend("bottomleft", c("t test", "Sign test"),fill=c("blue","red"))
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
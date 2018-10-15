plot(spbdry, axes=TRUE)
plot(sproads, add=TRUE, lty=2)
plot(spasthma, add=TRUE, pch=c(4,17)[(spasthma$Asthma == "case") + 1], cex=c(0.6, 0.75)[(spasthma$Asthma == "case") + 1])
plot(spsrc, pch=22, add=TRUE, cex=1.2, bg="grey60")



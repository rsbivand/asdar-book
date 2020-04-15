oopar <- par(mfrow=c(1,2), mar=c(3,3,1,1)+0.1)
plot(auck_shore)
legend("bottomleft", legend="a)", bty="n")
plot(auck_shore)
plot(islands_sp, add=TRUE, col="grey")
legend("bottomleft", legend="b)", bty="n")
par(oopar)



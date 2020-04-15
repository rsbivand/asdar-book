ylim <- range(c(nysarw_ll$llprof$ll, nysar_ll$llprof$ll), na.rm=TRUE)
plot(nysarw_ll$llprof$lambda, nysarw_ll$llprof$ll, type="l", xlab=expression(lambda), ylab="log likelihood", ylim=ylim, lwd=2)
abline(v=nysarw_ll$lambda)
abline(h=nysarw_ll$LL)
lines(nysar_ll$llprof$lambda, nysar_ll$llprof$ll, lty=2, lwd=2)
abline(v=nysar_ll$lambda, lty=2)
abline(h=nysar_ll$LL, lty=2)
legend("bottom", legend=c("weighted SAR", "SAR"), lty=c(1,2), lwd=2, bty="n")



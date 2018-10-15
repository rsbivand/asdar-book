plot(s, khcases-khcontrols, type="l", 
  ylab="D(s)", ylim=c(-.015, .015))#ylim=c(-11.5, 11.5))
lines(s, -1.96*sqrt(diag(khcov)), lty=2)
lines(s, +1.96*sqrt(diag(khcov)), lty=2)

envel<-apply(khcasesrel-khcontrolsrel, 1, function(X){quantile(X, c(.025, .975))})
lines(s, envel[1,], lty=3, lwd=2)
lines(s, envel[2,], lty=3, lwd=2)

legend("bottomleft", 
   legend=c("Actual value", "Approx. 95% C.I.", "Sim. 95% envelopes"),
   lty=1:3, lwd=c(1,1,2), bty="n")



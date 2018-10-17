plot(k, ylab="Intensity", main="")
points(x, rep(0, nx), pch=20)
for(i in 1:length(x))
	lines(density(x[i], bw=bw, kernel="biweight"), lty=2)

legend(x=14, y=0.6, legend=c("Intensity", "Kernel"), lty=c(1,2))



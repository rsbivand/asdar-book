def.par <- par(no.readonly = TRUE)
oopar <- par(mar=c(3,3,3,1)+0.1, mfrow=c(1,2))
plot(q5, pal=pal, main="Quantile", xlab="", ylab="")
plot(fj5, pal=pal, main="Fisher-Jenks", xlab="", ylab="")
par(oopar)
par(def.par)
cat("\n")



oopar <- par(mar=c(4,4,1,1)+0.1, mfrow=c(2,1))
plot(coordinates(transect_el1)[,1], transect_el1$band1, type="l", main="", xlab="", ylab="elevation, m", axes=FALSE)
abline(h=0)
box()
axis(2)
axis(1, at=axTicks(1), labels=parse(text=sp:::degreeLabelsEW(axTicks(1))))
plot(ecdf(transect_el1$band1), verticals= TRUE, do.p = FALSE, main="", xlab="elevation, m", ylab="")
par(oopar)



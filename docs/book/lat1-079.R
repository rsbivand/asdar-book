oopar <- par(mfrow=c(1,2))
msp <- moran.plot(NY8$Cases, listw=nb2listw(NY_nb, style="C"), quiet=TRUE)
title("Moran scatterplot")
infl <- apply(msp$is.inf, 1, any)
x <- NY8$Cases
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L", "H"), include.lowest=TRUE)
wx <- lag(nb2listw(NY_nb, style="C"), NY8$Cases)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)), labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4
plot(NY8, col=grey.colors(4, 0.95, 0.55, 2.2)[cols])
legend("topright", legend=c("None", "HL", "LH", "HH"), fill=grey.colors(4, 0.95, 0.55, 2.2), bty="n", cex=0.8, y.intersp=0.8)
title("Tracts with influence")
par(oopar)



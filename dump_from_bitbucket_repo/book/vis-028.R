def.par <- par(no.readonly = TRUE)
pin <- par("pin")
dxy <- apply(bbox(meuse),1,diff)
ratio <- dxy[1]/dxy[2]
par(pin = c(ratio * pin[2], pin[2]))
par(xaxs="i")
par(yaxs="i")
plot(meuse, pch = 1)
box()
par(def.par)
cat("\n")



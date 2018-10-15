library(lattice)
trellis.par.set(canonical.theme(color = FALSE))
library(RColorBrewer)
gry <- grey.colors(9, 0.95, 0.55, 2.2)
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(NY8, "lmresid", sp.layout=list(TCEpts), col.regions=colorRampPalette(gry)(7), at=seq(-2.5,5,length.out=6))



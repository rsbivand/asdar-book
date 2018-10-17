library(RColorBrewer)
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(NY8, "lmwresid", sp.layout=list(TCEpts), col.regions=colorRampPalette(gry)(7), at=seq(-2.5,5,length.out=6))



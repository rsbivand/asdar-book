library(RColorBrewer)
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(NY8, "lmpresid", sp.layout=list(TCEpts), col.regions=grey.colors(9, 0.95, 0.55, 2.2), at=seq(-3,5,length.out=9))



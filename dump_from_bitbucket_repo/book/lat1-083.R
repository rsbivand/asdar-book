library(lattice)
trellis.par.set(canonical.theme(color = FALSE))
library(RColorBrewer)
NY8$Standard <- lm1[,1]
NY8$"Constant_risk" <- I_CR
nms <- match(c("Standard", "Constant_risk"), names(NY8))
spplot(NY8, nms, at=c(-2.5,-1.4,-0.6,-0.2,0,0.2,0.6,4,7), col.regions=grey.colors(8, 0.95, 0.55, 2.2))



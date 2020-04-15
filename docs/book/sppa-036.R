grd <- GridTopology(cellcentre.offset=c(0.005,0.005), cellsize=c(0.01, 0.01),
  cells.dim=c(100, 100))
lambda<-exp(apply(coordinates(grd),1, function(X, alpha, beta)
	{
		loglambda(X, alpha, beta)
	}, alpha=optbeta$par[1], beta=optbeta$par[-1]	
 ))

parint<-SpatialGridDataFrame(grd, data=data.frame(intensity=lambda))

lyt<-list("sp.points", SpatialPoints(x), pch=19, col="black", cex=0.7)
print(spplot(parint, at=seq(0,1400,length.out=8),
 col.regions=colorRampPalette(gp)(7), sp.layout=list(lyt)))



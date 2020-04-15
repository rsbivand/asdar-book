cuts=c(0,20,50,200,500,2000)
grys <- grey.colors(7, 0.90, 0.50, 2.2)
print(spplot(meuse[1:4], main = "ppm", cuts=cuts, cex=.5, col.regions=grys),
      split=c(1,1,2,1),more=T)
meuse$lead.st = as.vector(scale(meuse$lead))
meuse$zinc.st = as.vector(scale(meuse$zinc))
meuse$copper.st = as.vector(scale(meuse$copper))
meuse$cadmium.st = as.vector(scale(meuse$cadmium))
cuts=c(-1.2,0,1,2,3,5)
print(spplot(meuse, c("lead.st", "zinc.st", "cadmium.st", "copper.st"),
	main = "standardised", cex = .5, cuts = cuts, col.regions=grys),
        split=c(2,1,2,1))
cat("\n")



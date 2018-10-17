oopar <- par(mar=c(1,1,1,1)+0.1)
library(RColorBrewer)
gcols <- grey.colors(15, 0.95, 0.55, 2.2)
image(sohoSG, "snowcost_broad", breaks=seq(0,750,50),
 col=gcols)
plot(buildings, col="white", add=TRUE)
plot(buildings, angle=45, density=10, col="grey70", add=TRUE)
symbols(coordinates(deaths), circles=4*sqrt(deaths$Num_Cases),
 inches=FALSE, add=TRUE, bg=c("grey75","grey50")[deaths$b_nearer+1])
source("legend_image.R") #from geoR
rect(528900, 180550, 529040, 180990, border=NA, col="white")
text(528970, 180950, "metres from\nBroad Street\npump", cex=0.6)
legend_image(c(528930, 528960), c(180600, 180900),
 sohoSG$snowcost_broad, vertical=TRUE, breaks=seq(0,750,50),
 col=gcols)
plot(nb_pump, add=TRUE, pch=8, cex=1.3, lwd=2)
plot(b_pump, add=TRUE, pch=4, cex=1.5, lwd=8, col="white")
plot(b_pump, add=TRUE, pch=4, cex=1.5, lwd=6)
rect(528900, 181330, 529140, 181380, border=NA, col="white")
legend(c(528910, 529100), c(181350, 181380),
 legend=c("Broad Street pump","other pumps"), pch=c(4,8), bty="n",
 cex=0.6, y.inter=0.7)
rect(528900, 181270, 529180, 181335, border=NA, col="white")
legend(c(528910, 529100), c(181275, 181325),
 legend=c("nearer Broad Street pump","nearer other pump"),
 fill=c("grey50","grey75"), bty="n", cex=0.6, y.inter=0.7)
box()
par(oopar)



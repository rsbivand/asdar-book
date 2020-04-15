library(RColorBrewer)
gp <- grey.colors(5, 0.9, 0.45, 2.2)


print(spplot(kernels, at=seq(0,2000,length.out=22),
 col.regions=colorRampPalette(gp)(21), 
names.attr=c(paste("Quartic bw=",bw, sep="", collapse=""),
"Quartic bw=0.05", "Quartic bw=0.1","Quartic bw=0.15", 
paste("Gaussian bw=",.5*bw, sep="", collapse=""),
 "Gaussian bw=0.025", "Gaussian bw=0.05","Gaussian bw=0.075") ) )






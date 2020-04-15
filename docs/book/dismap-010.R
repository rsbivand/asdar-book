library(RColorBrewer)

logSMR<-log(nc$SMR[nc$SMR>0])
nsteps <- 5
step<-(max(logSMR)-min(logSMR))/nsteps
brks<-exp(min(logSMR)+(0:nsteps)*step)
brks[1]<-0
cols <- grey.colors(nsteps, 0.95, 0.55, 2.2)
grps<-as.ordered(cut(nc$SMR, brks, include.lowest=TRUE))
plot(nc, col=cols[unclass(grps)], axes = FALSE)
box()
degAxis(1)
degAxis(2, at=c(34,35,36,37)) 
legend("bottomleft",legend=levels(grps), fill=cols, bty="n",cex=0.8,y.intersp=0.8) 



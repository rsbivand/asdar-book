library(RColorBrewer)
cl05 <- cl[cl$level == "0.05",]
xzx <- slot(slot(cl05, "lines")[[1]], "Lines")
cl05a <- SpatialLines(list(Lines(xzx, ID="0.05")))
lyt05 <- list("sp.lines", cl05a, lwd=2, lty=2, col="grey95")
lyt95 <- list("sp.lines", cl[cl$level == "0.95",], lwd=2, lty=1)
lytb <- list("sp.polygons", spbdry)
lytp <- list("sp.points", spsrc, cex=0.9, pch=4, col="grey95", lwd=3)
brks <- quantile(spkratio$kratio[spkratio$kratio>0], seq(0,1,1/10), na.rm=TRUE)
brks[1] <- 0
lbrks <- formatC(brks, 3, 6, "g", " ")
cols <- colorRampPalette(grey.colors(5, 0.95, 0.55, 2.2))(length(brks)-1)
colorkey<-list(labels=lbrks,
  at=(0:10)/10, height=.5)

print(spplot(spkratio, "kratio",
   col.regions=cols,
   do.log=TRUE, 
   colorkey=colorkey,
   at=c(0, brks[-c(1,11)], max(spkratio$kratio, na.rm=TRUE)),
   sp.layout=list(lyt05, lyt95, lytb, lytp) 
))



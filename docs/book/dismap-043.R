vbrks<-quantile(nc$BYMvmean, seq(0,1,1/5))
vbrks[6] <-vbrks[6]*1.01
colorkey$labels<-as.character(trunc(1000*vbrks)/1000)
colorkey$at<-vbrks[1]+vbrks[6]*0:5/5
print(spplot(nc, "BYMvmean", at=vbrks, axes=TRUE,
   col.regions=grey.colors(5, 0.95, 0.55, 2.2), colorkey=colorkey))



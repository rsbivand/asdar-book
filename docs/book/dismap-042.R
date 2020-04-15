ubrks<-quantile(nc$BYMumean, seq(0,1,1/5))
ubrks[6] <- ubrks[6]*1.01
colorkey$labels<-as.character(trunc(1000*ubrks)/1000)
colorkey$at<-ubrks[1]+ubrks[6]*0:5/5
print(spplot(nc, "BYMumean", at= ubrks, axes=TRUE, 
col.regions=grey.colors(5, 0.95, 0.55, 2.2), colorkey=colorkey))




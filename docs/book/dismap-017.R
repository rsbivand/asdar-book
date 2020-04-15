atcol<-(0:5)*max(nc$SMR)/5
colorkey<-list(labels=as.character(c(formatC(brks, format="f", dig=2))),
  at=atcol,  height=.5)

cols <- grey.colors(5, 0.95, 0.45, 2.2)

print(spplot(nc, c("SMR","EBPG", "EBLN", "EBMarshall"), col.regions=cols, 
  at=brks, axes = TRUE, colorkey=colorkey))




clusters<-get.knclusters(as(nc, "data.frame"), knresults)
i<-which.max(knresults$statistic)

nc$KNcluster<-""
nc$KNcluster[clusters[[i]]]<-"cluster"
nc$KNcluster[clusters[[i]][1]]<-"centre"
nc$KNcluster<-as.factor(nc$KNcluster)

print(spplot(nc, "KNcluster", main="Kulldorff's method",
  xlab="Easting", ylab="Northing", col.regions=c(gray(1), gray(.5), gray(.8))))





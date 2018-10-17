oopar <- par(mfrow=c(1,2), mar=c(5,3,1,1)+0.1)
b_wid <- table(deaths$b_nearer)
boxplot(snowcost_broad ~ b_nearer, deaths, width=b_wid, ylim=c(0,450), ylab="distance", xlab="Broad Street", col=grey.colors(1, 0.8, 0.8, 2.2))
boxplot(snowcost_not_broad ~ b_nearer, deaths, width=b_wid, ylim=c(0,450), xlab="Other pump", col=grey.colors(1, 0.8, 0.8, 2.2))
par(oopar)



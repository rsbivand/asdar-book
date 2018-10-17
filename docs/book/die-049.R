oopar <- par(mar=c(4,7,2,2)+0.1, las=1)
tg <- table(spear$geology)
boxplot(elevation.dem ~ geology, spear, width=tg, horizontal=TRUE)
par(oopar)



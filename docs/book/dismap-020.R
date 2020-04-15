oopar <- par(mar=c(3,7,2,1)+0.1)
boxplot(as(nc, "data.frame")[,c("SMR", "EBPG", "EBLN", "EBMarshall", "EBMrshloc")], cex.lab=.5, las=1, horizontal=TRUE)
par(oopar)



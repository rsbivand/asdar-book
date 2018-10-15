oopar <- par(mar=c(1,1,1,1)+0.1)
plot(manitoulin_sp, pbg="grey75", col="grey95")
text(t(sapply(slot(slot(manitoulin_sp, "polygons")[[1]], "Polygons"), function(x) slot(x, "labpt")))[-c(1,2),], label=high$polydata$level[-c(1,2)], col="black", font=2)
par(oopar)



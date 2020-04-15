library(lattice)
print(levelplot(z~x+y|name, spmap.to.lev(zn[c("direct", "log")]), asp = "iso",
        cuts=4, col.regions=grey.colors(5, 0.90, 0.50, 2.2)),
	split = c(1,1,1,2), more = TRUE)
print(spplot(zn[c("direct", "log")], cuts=4,
        col.regions=grey.colors(5, 0.90, 0.50, 2.2)), split = c(1,2,1,2))
cat("\n")



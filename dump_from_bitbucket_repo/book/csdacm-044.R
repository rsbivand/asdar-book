library(lattice)
grys <- grey.colors(11, 0.95, 0.55, 2.2)
print(spplot(meuse.grid["dist"], cuts=10, col.regions=grys, sp.layout = list("sp.points", HexPts, col = 1)),
	split = c(1, 1, 2, 1), more = TRUE)
print(spplot(HexPolsDf["dist"], cuts=10, col.regions=grys),
	split = c(2, 1, 2, 1), more = FALSE)



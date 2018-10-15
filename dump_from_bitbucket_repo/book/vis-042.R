grys <- grey.colors(7, 0.90, 0.50, 2.2)
print(spplot(zn["log"], sp.layout = meuse.layout, cuts=6, col.regions=grys))
cat("\n")



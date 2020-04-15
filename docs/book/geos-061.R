x$lz.uk <- lz.uk$var1.pred
x$lz.ok <- lz.ok$var1.pred
print(spplot(x, c("log.zinc.pred", "lz.ok", "lz.uk"),
	names.attr = c("collocated", "ordinary", "universal"),
	cuts=7, col.regions=grey.colors(8, 0.55, 0.95, 2.2)
))



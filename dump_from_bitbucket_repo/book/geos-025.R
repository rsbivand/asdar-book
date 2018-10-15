sel <-
structure(list(x = c(145.291968730077, 266.107479142605, 320.156523274526, 
339.232656497557, 323.335878811698, 212.058435010685, 135.753902118561, 
46.7319470777507, 78.5255024494688, 142.112613192905), y = c(574649.690841889, 
581256.265954825, 627502.29174538, 822396.257577002, 1053626.38652977, 
1278249.94036961, 1255126.92747433, 792666.669568789, 634108.866858316, 
577952.978398357)), .Names = c("x", "y"))
v <- variogram(zinc ~ 1, meuse, cloud = TRUE)
v$gamma <- v$gamma/1e6
sel$y <- sel$y/1e6
p1 <- xyplot(gamma~dist, v, 
	panel = function(x, y, ...) {
		panel.xyplot(x, y, ...)
		llines(sel$x, sel$y)
	},
	pch=3, cex = .5, asp = 1, ylab = "gamma (x 1e6)")
x <-
structure(list(head = c(40, 40, 40, 54, 55, 54, 47, 80, 55, 55, 
54, 53, 54, 55, 59, 59), tail = c(41, 42, 43, 57, 57, 58, 59, 
99, 121, 122, 123, 125, 125, 125, 125, 132)), .Names = c("head", 
"tail"), row.names = as.integer(c(NA, 16)), class = c("pointPairs", 
"data.frame"))
p2 = plot(x, meuse, scales=list(draw=F), col.line = 1)
print(p1, split = c(1,1,2,1), more = TRUE)
print(p2, split = c(2,1,2,1))



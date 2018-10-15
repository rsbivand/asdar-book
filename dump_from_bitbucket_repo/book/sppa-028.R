plot(mserw$h, mserw$mse, xlab="Bandwidth", ylab="MSE", type="l", ylim=c(-2,50))
i<-which.min(mserw$mse)
points(mserw$h[i], mserw$mse[i])



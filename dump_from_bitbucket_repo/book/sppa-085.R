oopar <- par(mfrow=c(1,2))
plot(s, kihnocov-kihmeannocov, type="l", 
   ylim= c(-0.06,  0.22),
   xlab="s", ylab= expression(hat(k)[I][","][hat(lambda)]-"E[s]"),
    main ="No covariates" )

envnocov<-apply(kinhomrelnocov, 1, function(X){quantile(X, c(.025, .975))})
lines(s, envnocov[1,]-kihmeannocov, lty=2)
lines(s, envnocov[2,]-kihmeannocov, lty=2)
plot(s, kih-kihmean, type="l", ylim=c(-0.06,  0.22), #c(-2e-4, 2e-4),
   xlab="s", ylab= expression(hat(k)[I][","][hat(lambda)]-"E[s]"),
   main ="Adjusting for Hay Fever"  )

env<-apply(kinhomrel, 1, function(X){quantile(X, c(.025, .975))})
lines(s, env[1,]-kihmean, lty=2)
lines(s, env[2,]-kihmean, lty=2)
par(oopar)



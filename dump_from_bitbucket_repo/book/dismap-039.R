library(coda)

plot(ncoutput[,c("deviance", "alpha", "beta", "theta[94]")])

     save(file="BYM.RData", list=c("d", "inits", "MCMCres", "ncoutput") )




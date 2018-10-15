ats <- seq(0,max(spkratio$prob),length.out=11)
cols <- colorRampPalette(grey.colors(5, 0.9, 0.5, 2.2))(length(ats)-1)
print(spplot(spkratio, "prob", col.regions=cols, at=ats, sp.layout=list(lytb, lytp)))



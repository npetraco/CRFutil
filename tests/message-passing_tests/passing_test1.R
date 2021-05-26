library(CRFutil)
library(Rgraphviz)

# Model:
grf <- ~A:B
plot(ug(grf))
dev.off()

# Convert to pw factor graph
pwfg <- mrf2pwfg2(grf, plotQ=T)
dev.off()
plot(pwfg, nodeAttrs=makeNodeAttrs(pwfg, fontsize=30))


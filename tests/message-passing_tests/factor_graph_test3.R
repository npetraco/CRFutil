library(CRFutil)

# Model:
#grf <- ~A:B
grf <- ~A:C + B:D + E:C + E:B
#grf <- ~1:2 # Node names must be numbers for now.... FIX
#grf <- ~1:3 + 2:3 + 3:4 + 3:5 # Node names must be numbers for now.... FIX
#grf <- ~1:2 + 2:4 + 3:4 + 4:5 + 4:6 # Node names must be numbers for now.... FIX
#grf <- ~1:2 + 1:4 + 2:5 + 2:6 + 3:7 + 4:7 + 4:8 + 6:9 + 7:10 + 11:10 + 10:14 + 10:13 + 10:12 + 12:17 + 12:16 + 12:15
grf <-  ~A:B + A:D + B:E + B:FF + C:G + D:G + D:H + FF:II + G:J + K:J + J:N + J:M + J:L + L:Q + L:P + L:O
plot(ug(grf))
dev.off()

pwfg <- mrf2pwfg(grf, plotQ=T)
plot(pwfg)
dev.off()

pwfg2 <- mrf2pwfg2(grf, plotQ=T) # This version can have character names for the nodes
library(Rgraphviz)
plot(pwfg2, nodeAttrs=makeNodeAttrs(pwfg2, fontsize=50))

tv<-c("A","B","C")
tv[c(3,2,1,3,1,1,3,1)]

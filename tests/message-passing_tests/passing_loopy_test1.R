library(CRFutil)
library(Rgraphviz)


# Model:
#grf <- ~1:2
#grf <- ~A:B
grf <- ~1:2 + 2:3 + 3:1
#grf <- ~A:B + B:C + C:A
#grf <- ~A:B + B:C + C:A + B:D
#grf <- ~A:B + B:C + C:A + B:FF + FF:G + G:B
#grf <- ~A:B + B:C + C:A + B:D + D:FF + FF:G + G:D
#grf <- ~A:B + B:C + C:A + B:E + E:D + D:FF + FF:G + G:D
#grf <- ~A:B + B:C + C:A + B:E + E:D + D:FF + FF:G + G:D + E:H

dev.off()
plot(ug(grf))
dev.off()


# Convert MRF to a pair-wise factor graph
pwfg <- mrf2pwfg(grf, plotQ=T)
dev.off()
plot(pwfg, nodeAttrs=makeNodeAttrs(pwfg, fontsize=30))

msg.cont.info <- init.loopy.message.storage(pwfg)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

msg.bxs
msg.sch

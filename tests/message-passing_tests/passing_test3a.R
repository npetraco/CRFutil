library(CRFutil)
library(Rgraphviz)

# NOTE!!!!! mrf2pwfg2 has replaced mrf2pwfg from this test script on!!!!!!!!


# Model:
grf <- ~1:2 + 1:4 + 2:5 + 2:6 + 3:7 + 4:7 + 4:8 + 6:9 + 7:10 + 11:10 + 10:14 + 10:13 + 10:12 + 12:17 + 12:16 + 12:15
#grf <-  ~A:B + A:D + B:E + B:FF + C:G + D:G + D:H + FF:II + G:J + K:J + J:N + J:M + J:L + L:Q + L:P + L:O
plot(ug(grf))
dev.off()

# First Instantiate a model with some random values for the params in the potentials:
adj      <- ug(grf, result="matrix")
n.states <- 2
km       <- sim.field.random(adjacentcy.matrix = adj, num.states = n.states, num.sims = 0)[[1]]

# Dress potentials for message passing:
gr.pots  <- make.gRbase.potentials(crf = km, node.names = rownames(adj), state.nmes = c("up", "dn"))
ep       <- gr.pots$edge.potentials # Edge potentials in gRbase format
np       <- gr.pots$node.potentials # Node potentials in gRbase format
all.pots <- union(np, ep)           # Put all node and edge potentials into one list so we don't have to search though two lists for them

nd.nms          <- nodes(ug(grf))                                               # node names
edg.nms         <- cbind(nd.nms[km$edges[,1]], nd.nms[km$edges[,2]])            # edge names
nd.pot.nms      <- paste0("f", nd.nms)                                          # node potential names
edg.pot.nms     <- paste0("f", nd.nms[km$edges[,1]], "-", nd.nms[km$edges[,2]]) # edge potential names
names(all.pots) <- c(nd.pot.nms, edg.pot.nms) # Put the same names on the potentials as will appear in the schedule
all.pots        # Potentials dressed up in gRbase format and ready for message passing

# Convert MRF to a pair-wise factor graph
pwfg <- mrf2pwfg(grf, plotQ=F)
dev.off()
plot(pwfg, nodeAttrs=makeNodeAttrs(pwfg, fontsize=30))

# Initialize a storage list to hold messages and get message passing schedule:
root.pths     <- get.root.paths(pwfg, root.node = "1", serial.schedsQ = T)
msg.cont.info <- init.message.storage(root.pths)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

# Pass messages according to schedule:
for(i in 1:nrow(msg.sch)){
  msg.num <- i
  print(msg.sch[msg.num,]) # message info
  msg.nd.nmes <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message

  msg.info <- pass.message(start.node = msg.nd.nmes[1],
                           end.node   = msg.nd.nmes[2],
                           factor.graph = pwfg,
                           pots.list  = all.pots,
                           mailboxes.list = msg.bxs,
                           printQ = F)

  msg.bxs[[msg.info$mailbox.idx]] <- msg.info$message

}
msg.bxs

# Get exact marginals to check against:
all.marginals <- infer.exact(km)
all.marginals

# Node marginals
nde.nms <- nodes(ug(grf))
for(i in 1:length(nde.nms)){
  ndm <- node.marginal(nde.nms[i], pwfg, msg.bxs)
  print(paste0("bel(",nde.nms[i],"):"))
  print(ndm)
  print(round(ndm - all.marginals$node.bel[i,], 3)) # Diffs = about 0??
  print("==================================")
}


# Edge marginals
for(i in 1:nrow(km$edges)){
  lnd <- nd.nms[km$edges[i,1]]
  rnd <- nd.nms[km$edges[i,2]]
  print(paste0("bel(",lnd,",",rnd,"):"))
  edgm <- edge.marginal(v.start.node = lnd, v.end.node = rnd,
                        factor.graph = pwfg, pots.list = all.pots, mailbox.list = msg.bxs)
  print(edgm)
  print(round(edgm - all.marginals$edge.bel[[i]], 3)) # Diffs = about 0??
  print("==================================")
}

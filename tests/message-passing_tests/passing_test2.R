library(CRFutil)
library(Rgraphviz)

# Model:
grf <- ~A:B
plot(ug(grf))
dev.off()

# Convert to pw factor graph
pwfg <- mrf2pwfg2(grf, plotQ=T) # Problem !!!!!!!!
dev.off()
plot(pwfg, nodeAttrs=makeNodeAttrs(pwfg, fontsize=30))

# Message pass test:

# First Instantiate a model with some random values for the params in the potentials:
adj      <- ug(grf, result="matrix") # NOTE: nodes are ordered according to how they show up in adj!!
n.states <- 2
knm      <- sim.field.random(adjacentcy.matrix = adj, num.states = n.states, num.sims = 0)[[1]]

knm$par # "true" theta
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)
gr.pots <- make.gRbase.potentials(crf = knm, node.names = rownames(adj), state.nmes = c("up", "dn"))

ep <- gr.pots$edge.potentials # Edge potentials in gRbase format
np <- gr.pots$node.potentials # Node potentials in gRbase format



# Initialize a storage list to hold messages and get message passing schedule:
root.pths     <- get.root.paths(pwfg, root.node = "A", serial.schedsQ = T)
msg.cont.info <- init.message.storage(root.pths)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

msg.bxs
msg.sch

msg.num <- 3
msg.sch[msg.num,]
msg.sp  <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message
msgs.in <- get.incoming.messages(start.node = msg.sp[1], end.node = msg.sp[2], factorgraph = pwfg, message.list = msg.bxs)
msgs.in


# Picking out the leaf (node) potentials:
# Name check incoming message for length and a name. Length 1 and no name indicates no incoming neighbors
rownames(adj)
np
ep[[2]]
msg.sp
make.f2v.msg(in.v.msgs.list = NULL, f.msg = )

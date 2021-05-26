library(CRFutil)

# Model:
#grf <- ~1:2# Node names must be numbers for now.... FIX
#grf <- ~1:3 + 2:3 + 3:4 + 3:5 # Node names must be numbers for now.... FIX
grf <- ~1:2 + 2:4 + 3:4 + 4:5 + 4:6 # Node names must be numbers for now.... FIX

# Convert to pw factor graph
pwfg <- mrf2pwfg(grf, plotQ=T)

# Initialize a storage list to hold messages and get message passing schedule:
root.pths     <- get.root.paths(pwfg, root.node = 4, serial.schedsQ = T)
msg.cont.info <- init.message.storage(root.pths)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

msg.bxs
msg.sch

# Grab messages of requisite neighboring nodes needed to construct the outgoing message
# This task needs msg.cont.info and the factor graph
msg.num      <- 16
msg.sch[msg.num,]
msg.sp       <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message
get.incoming.messages(start.node = msg.sp[1], end.node = msg.sp[2], factorgraph = pwfg, message.list = msg.bxs)

# Message pass test:
# First Instantiate a model with some random values for the params in the potentials:
adj      <- ug(grf, result="matrix")
n.states <- 2
knm      <- sim.field.random(adjacentcy.matrix = adj, num.states = n.states, num.sims = 0)[[1]]

knm$par # "true" theta
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)

rownames(adj) # NOTE Node names probably aren't in consecutive order!!!!!!!!!

gr.pots <- make.gRbase.potentials(crf = knm, node.names = rownames(adj), state.nmes = c("up", "dn"))
ep <- gr.pots$edge.potentials # Edge potentials in gRbase format
np <- gr.pots$node.potentials # Node potentials in gRbase format
np
ep[[2]]
knm$edge.pot[[2]]

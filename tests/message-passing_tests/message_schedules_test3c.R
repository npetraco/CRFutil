library(CRFutil)

# Model:
#grf <- ~1:2 # Node names must be numbers for now.... FIX
#grf <- ~1:3 + 2:3 + 3:4 + 3:5 # Node names must be numbers for now.... FIX
#grf <- ~1:2 + 2:4 + 3:4 + 4:5 + 4:6 # Node names must be numbers for now.... FIX
#grf <- ~1:2 + 1:4 + 2:5 + 2:6 + 3:7 + 4:7 + 4:8 + 6:9 + 7:10 + 11:10 + 10:14 + 10:13 + 10:12 + 12:17 + 12:16 + 12:15
grf <-  ~A:B + A:D + B:E + B:FF + C:G + D:G + D:H + FF:II + G:J + K:J + J:N + J:M + J:L + L:Q + L:P + L:O
plot(ug(grf))

# Convert to pw factor graph
pwfg <- mrf2pwfg2(grf, plotQ=T)
library(Rgraphviz)
dev.off()
plot(pwfg, nodeAttrs=makeNodeAttrs(pwfg2, fontsize=30))

# Initialize a storage list to hold messages and get message passing schedule:
root.pths     <- get.root.paths(pwfg, root.node = "J", serial.schedsQ = T)
msg.cont.info <- init.message.storage(root.pths)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

msg.bxs
msg.sch

# Grab messages of requisite neighboring nodes needed to construct the outgoing message
# This task needs msg.cont.info and the factor graph
for(i in 1:nrow(msg.sch)) {

  msg.num <- i
  msg.sp  <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message
  slj     <- get.incoming.messages(start.node = msg.sp[1], end.node = msg.sp[2], factorgraph = pwfg, message.list = msg.bxs)

  print(names(slj))
  print(length(names(slj)))
  print("-------------------")

}


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

msg.sch
msg.num <- 1
msg.sch[msg.num,]
msg.sp  <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message
msgs.in <- get.incoming.messages(start.node = msg.sp[1], end.node = msg.sp[2], factorgraph = pwfg, message.list = msg.bxs)
msgs.in

# Put names on gr node potentials
# Put names on gr edge potentials??
# Name check incoming message for length and a name. Length 1 and no name indicates no incoming neighbors
rownames(adj)
np
ep[[2]]
msg.sp
make.f2v.msg(in.v.msgs.list = NULL, f.msg = )



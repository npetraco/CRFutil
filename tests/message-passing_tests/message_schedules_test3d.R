library(CRFutil)
library(Rgraphviz)

# Model:
#grf <- ~1:2 # Node names must be numbers for now.... FIX
#grf <- ~1:3 + 2:3 + 3:4 + 3:5 # Node names must be numbers for now.... FIX
#grf <- ~1:2 + 2:4 + 3:4 + 4:5 + 4:6 # Node names must be numbers for now.... FIX
#grf <- ~1:2 + 1:4 + 2:5 + 2:6 + 3:7 + 4:7 + 4:8 + 6:9 + 7:10 + 11:10 + 10:14 + 10:13 + 10:12 + 12:17 + 12:16 + 12:15
grf <-  ~A:B + A:D + B:E + B:FF + C:G + D:G + D:H + FF:II + G:J + K:J + J:N + J:M + J:L + L:Q + L:P + L:O
plot(ug(grf))
dev.off()

# Convert to pw factor graph
pwfg <- mrf2pwfg2(grf, plotQ=T)
dev.off()
plot(pwfg, nodeAttrs=makeNodeAttrs(pwfg, fontsize=30))

# Initialize a storage list to hold messages and get message passing schedule:
root.pths     <- get.root.paths(pwfg, root.node = "J", serial.schedsQ = T)
msg.cont.info <- init.message.storage(root.pths)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

msg.bxs
msg.sch


# Message pass test:
# First Instantiate a model with some random values for the params in the potentials:
adj      <- ug(grf, result="matrix") # NOTE: nodes are ordered according to how they show up in adj!!
n.states <- 2
knm      <- sim.field.random(adjacentcy.matrix = adj, num.states = n.states, num.sims = 0)[[1]]

knm$par # "true" theta
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)

rownames(adj) # NOTE Node names probably aren't in consecutive order!!!!!!!!!
nodes(ug(grf))

gr.pots <- make.gRbase.potentials(crf = knm, node.names = rownames(adj), state.nmes = c("up", "dn"))

# NEED TO CHECK THAT EDGE POTENTIAL NAMES ARE PROPERLY ALIGNED WITH MESSAGE NAMES??
#TEST:
#Get same result with:
# a b  vs a c  ?? ie fA-B VS fB-A?? IE, do we have to account if edge potentials are asymetric??
# c d     b d
knm$edges
knm$node.par

nd.nms  <- nodes(ug(grf))
edg.nms <- cbind(nd.nms[knm$edges[,1]], nd.nms[knm$edges[,2]])

ep <- gr.pots$edge.potentials # Edge potentials in gRbase format
np <- gr.pots$node.potentials # Node potentials in gRbase format


# Put names on gr node potentials??  **Don't think we need too as long as we use the same adj mat for the graph and to crf object.
# Put names on gr edge potentials??  ** Same

sapply(1:length(np), function(xx){names(dimnames(np[[xx]]))})
nd.nms # Same??

# col 1 == col 3, col 2 == col 4??
cbind(t(sapply(1:length(ep), function(xx){names(dimnames(ep[[xx]]))})), edg.nms)

# If so, we can use knm$edges row indices to retrieve the required edge potential from ep

# Also, check eg:
ep[[2]] - knm$edge.pot[[2]]
lapply(1:length(ep), function(xx){ep[[xx]] - knm$edge.pot[[xx]]})

# Potential names. Node name order within the edge potential name should be the same as in msg.sch output by init.message.storage
# and (more fundamentally) mrf2pwfg2. In mrf2pwfg2 we use node.names <- und.gph@nodes which sets the order downstream
# Put all node and edge potentials into one list so we dont have to search though two lists for them
all.pots <- union(np, ep)
length(np)
length(ep)
length(all.pots)

# Put the same names on the potentials as will appear in the schedule:
nd.nms  <- nodes(ug(grf))                                      # node names
edg.nms <- cbind(nd.nms[knm$edges[,1]], nd.nms[knm$edges[,2]])   # edge names
names(all.pots) <- c(paste0("f", nd.nms), paste0("f", nd.nms[knm$edges[,1]], ".", nd.nms[knm$edges[,2]]))
all.pots




msg.sch
msg.num <- 1
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



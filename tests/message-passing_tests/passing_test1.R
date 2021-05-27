library(CRFutil)
library(Rgraphviz)

# Message pass test

# Model:
grf <- ~A:B
adj <- ug(grf, result="matrix")
f0  <- function(y){ as.numeric(c((y==1),(y==2)))}

n.states <- 2
km       <- make.crf(adj, n.states)

# Node pots:
PsiA <- c(3,1)
PsiB <- c(3,1)

# Edge pots:
PsiAB <-
  rbind(
    c(3, 6.1),
    c(6.1, 3)
  )

km$node.pot[1,]  <- PsiA
km$node.pot[2,]  <- PsiB
km$edge.pot[[1]] <- PsiAB


# theta corresponding to the potentials above. Needed to dress the potentials with requise labeling:
#km$par <- make.par.from.all.potentials(km)$par # If there is no par in crf obj
#out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T) # If there are no potentials in crf object

#XXXXXXXXXXXXXXXXXXXXXX
# Below function: send in crf object and extract/prepare potentials for message passing routines:


gr.pots <- make.gRbase.potentials(crf = km, node.names = rownames(adj), state.nmes = c("up", "dn"))

ep <- gr.pots$edge.potentials # Edge potentials in gRbase format
np <- gr.pots$node.potentials # Node potentials in gRbase format

# Put all node and edge potentials into one list so we dont have to search though two lists for them
all.pots <- union(np, ep)

# Put the same names on the potentials as will appear in the schedule:
nd.nms  <- nodes(ug(grf))                                      # node names
edg.nms <- cbind(nd.nms[km$edges[,1]], nd.nms[km$edges[,2]])   # edge names
names(all.pots) <- c(paste0("f", nd.nms), paste0("f", nd.nms[km$edges[,1]], ".", nd.nms[km$edges[,2]]))
all.pots

#XXXXXXXXXXXXXXXXXXXXXX
# Function below: send in graph eq to create a corresponding pair-wise factor graph and a message
# passing schedule. Graph should have no loops!

# Convert MRF to a pair-wise factor graph
pwfg <- mrf2pwfg2(grf, plotQ=T)
dev.off()
plot(pwfg, nodeAttrs=makeNodeAttrs(pwfg, fontsize=30))

# Initialize a storage list to hold messages and get message passing schedule:
root.pths     <- get.root.paths(pwfg, root.node = "A", serial.schedsQ = T)
msg.cont.info <- init.message.storage(root.pths)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

msg.bxs # message storage list (mailboxes)
msg.sch # schedule


#XXXXXXXXXXXXXXXXXXXXXXXXXXX
# Pass messages according to schedule:
msg.num <- 1
msg.sch[msg.num,]
msg.typ <- message.type(msg.sch[msg.num,4])
msg.sp  <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message
msgs.in <- get.incoming.messages(start.node = msg.sp[1], end.node = msg.sp[2], factorgraph = pwfg, message.list = msg.bxs)
msgs.in

msg.sp
msg.typ
if(msg.typ == "f2v"){
  fpot <- all.pots[[which(names(all.pots) == msg.sp[1])]]
}
make.f2v.msg(in.v.msgs.list = NULL, f.msg = fpot)
make.f2v.msg(in.v.msgs.list = msgs.in, f.msg = fpot)


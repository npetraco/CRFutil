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
names(all.pots) <- c(paste0("f", nd.nms), paste0("f", nd.nms[km$edges[,1]], "-", nd.nms[km$edges[,2]]))
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
msg.num <- 8
msg.sch[msg.num,] # message info

# Get message name and the mailbox number it goes in:
msg.nme <- msg.sch[msg.num,4]
msg.idx <- which(names(msg.bxs) == msg.nme)
msg.nme
msg.idx

# Get incoming messages:
msg.sp  <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message
msgs.in <- get.incoming.messages(start.node = msg.sp[1], end.node = msg.sp[2], factorgraph = pwfg, message.list = msg.bxs)
msgs.in

# Determine is f2v of v2f type message
msg.typ <- message.type(msg.sch[msg.num,4])
msg.typ

msg.sp
msg.typ
msg.nme
msg.idx

if(msg.typ == "f2v"){
  fpot <- all.pots[[which(names(all.pots) == msg.sp[1])]]
  msg.bxs[[msg.idx]] <- make.f2v.msg(in.v.msgs.list = msgs.in, f.msg = fpot, out.v.nme = msg.sp[2])
} else if(msg.typ == "v2f"){
  msg.bxs[[msg.idx]] <- make.v2f.msg(in.f.msgs.list = msgs.in)
} else {
  stop("Invalid message pass requested!")
}
msg.bxs

# Checks:
msg.num <- 1      # message pass 1
msg.sch[msg.num,] # message info

# Get message name and the mailbox number it goes in:
msg.nme <- msg.sch[msg.num,4]
msg.idx <- which(names(msg.bxs) == msg.nme)
msg.nme
msg.idx

# Should be \Psi_A
names(all.pots)

msg.bxs[[msg.idx]]
all.pots[[1]]


#-----------------------------------------
msg.num <- 2      # message pass 1
msg.sch[msg.num,] # message info

# Get message name and the mailbox number it goes in:
msg.nme <- msg.sch[msg.num,4]
msg.idx <- which(names(msg.bxs) == msg.nme)
msg.nme
msg.idx

# Should be \Psi_B
names(all.pots)

msg.bxs[[msg.idx]]
all.pots[[2]]

#-----------------------------------------
msg.num <- 3      # message pass 2
msg.sch[msg.num,] # message info

# Get message name and the mailbox number it goes in:
msg.nme <- msg.sch[msg.num,4]
msg.idx <- which(names(msg.bxs) == msg.nme)
msg.nme
msg.idx

# Should be \Psi_B
names(all.pots)

msg.bxs[[msg.idx]]
all.pots[[2]]


#-----------------------------------------
msg.num <- 4      # message pass 3
msg.sch[msg.num,] # message info

# Get message name and the mailbox number it goes in:
msg.nme <- msg.sch[msg.num,4]
msg.idx <- which(names(msg.bxs) == msg.nme)
msg.nme
msg.idx

# Should be \sum_{\sim A} \Big( \Psi_{\text{AB}} \times \Psi_{B} \Big)
names(all.pots)
tabMarg(tabProd(all.pots[[3]], all.pots[[2]]), marg = "A")

msg.bxs[[msg.idx]]


#-----------------------------------------
msg.num <- 5      # message pass 4
msg.sch[msg.num,] # message info

# Get message name and the mailbox number it goes in:
msg.nme <- msg.sch[msg.num,4]
msg.idx <- which(names(msg.bxs) == msg.nme)
msg.nme
msg.idx

# Should be \sum_{\sim A} \Big( \Psi_{\text{AB}} \times \Psi_{B} \Big)
names(all.pots)
tabMarg(tabProd(all.pots[[3]], all.pots[[2]]), marg = "A")

msg.bxs[[msg.idx]]


#-----------------------------------------
msg.num <- 6      # message pass 4
msg.sch[msg.num,] # message info

# Get message name and the mailbox number it goes in:
msg.nme <- msg.sch[msg.num,4]
msg.idx <- which(names(msg.bxs) == msg.nme)
msg.nme
msg.idx

# Should be \Psi_{A}
names(all.pots)
all.pots[[1]]
msg.bxs[[msg.idx]]



#-----------------------------------------
msg.num <- 7      # message pass 5
msg.sch[msg.num,] # message info

# Get message name and the mailbox number it goes in:
msg.nme <- msg.sch[msg.num,4]
msg.idx <- which(names(msg.bxs) == msg.nme)
msg.nme
msg.idx

# Should be \sum_{\sim B} \Big( \Psi_{\text{AB}} \times \Psi_A \Big)
names(all.pots)
tabMarg(tabProd(all.pots[[3]], all.pots[[1]]), marg = "B")

msg.bxs[[msg.idx]]


#-----------------------------------------
msg.num <- 8      # message pass 6
msg.sch[msg.num,] # message info

# Get message name and the mailbox number it goes in:
msg.nme <- msg.sch[msg.num,4]
msg.idx <- which(names(msg.bxs) == msg.nme)
msg.nme
msg.idx

# Should be \sum_{\sim B} \Big( \Psi_{\text{AB}} \times \Psi_A \Big)
names(all.pots)
tabMarg(tabProd(all.pots[[3]], all.pots[[1]]), marg = "B")

msg.bxs[[msg.idx]]



library(CRFutil)
library(Rgraphviz)

# Message pass test

# Model:
grf <- ~A:B + B:C + C:D + C:E
adj <- ug(grf, result="matrix")
f0  <- function(y){ as.numeric(c((y==1),(y==2)))}
plot(ug(grf))
dev.off()

n.states <- 2
km       <- make.crf(adj, n.states)

# Node pots:
PsiA <- c(2,1)
PsiB <- c(3.5,1)
PsiC <- c(3,1)
PsiD <- c(2.5,1)
PsiE <- c(0.5,1)

# Edge pots:
PsiAB <-
  rbind(
    c(3, 6.1),
    c(6.1, 3)
  )
PsiBC <-
  rbind(
    c(6.1, 3),
    c(3, 6.1)
  )
PsiCD <-
  rbind(
    c(4, 6.2),
    c(6.2, 4)
  )
PsiCE <-
  rbind(
    c(3.5, 5.1),
    c(5.1, 3.5)
  )

km$node.pot[1,]  <- PsiA
km$node.pot[2,]  <- PsiB
km$node.pot[3,]  <- PsiC
km$node.pot[4,]  <- PsiD
km$node.pot[5,]  <- PsiE
km$edge.pot[[1]] <- PsiAB
km$edge.pot[[2]] <- PsiBC
km$edge.pot[[3]] <- PsiCD
km$edge.pot[[4]] <- PsiCE

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
pwfg <- mrf2pwfg2(grf, plotQ=F)
dev.off()
plot(pwfg, nodeAttrs=makeNodeAttrs(pwfg, fontsize=30))

# Initialize a storage list to hold messages and get message passing schedule:
root.pths     <- get.root.paths(pwfg, root.node = "D", serial.schedsQ = T)
msg.cont.info <- init.message.storage(root.pths)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

# Pass messages according to schedule:
# msg.num <- 3
# print(msg.sch[msg.num,]) # message info
# msg.nd.nmes <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message
#
# pass.message(start.node = msg.nd.nmes[1],
#              end.node   = msg.nd.nmes[2],
#              factor.graph = pwfg,
#              pots.list  = all.pots,
#              mailboxes.list = msg.bxs, printQ = F)

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


all.marginals <- infer.exact(km)

# node marginals ????
nde.nms <- nodes(ug(grf))
nde.idx <- 5
nde.nms[nde.idx]
node.marginal(nde.nms[nde.idx], pwfg, msg.bxs)
all.marginals$node.bel[nde.idx,] # Compare: same??
lapply(1:length(nde.nms), function(xx){node.marginal(nde.nms[xx], pwfg, msg.bxs)})

# edge marginals ??
edge.marginal(v.start.node = "E", v.end.node = "C", factor.graph = pwfg, pots.list = all.pots, mailbox.list = msg.bxs)
all.marginals$edge.bel[[4]]

for(i in 1:nrow(km$edges)){
  lnd <- nd.nms[km$edges[i,1]]
  rnd <- nd.nms[km$edges[i,2]]
  print(paste0("bel(",lnd,",",rnd,"):"))
  edgm <- edge.marginal(v.start.node = lnd, v.end.node = rnd,
                        factor.graph = pwfg, pots.list = all.pots, mailbox.list = msg.bxs)
  print(edgm)
  print(all.marginals$edge.bel[[i]])
}

# make identity message for loopy initalization
# Bethe
